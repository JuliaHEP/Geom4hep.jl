#---Meshes with normals for drawing------------------------------------------------------------------

"""
    shape_to_mesh(shape; facets=64, transformation=one(Transformation3D{T}), crease_angle=π/6)

Surface mesh of `shape` (a `GeometryBasics.Mesh` of triangles) with per-vertex normals, for smooth
shading of curved surfaces and sharp shading of edges.

The mesh of [`mesh`](@ref) is processed as follows:
- duplicated vertices are merged, and the faces are split into triangles;
- the triangles are oriented outwards, using the [`normal`](@ref) of the shape at their vertices
  (and, where the vertices are not on the surface, their neighbouring triangles);
- the normal of a vertex is the average of the normals of the triangles around it that make an
  angle smaller than `crease_angle` with each other. A vertex on a hard edge is duplicated, one
  copy for each side of the edge. When the normal of the shape at a vertex agrees with this
  average, it is used instead, which is exact on curved surfaces.

`facets` is the number of segments of the curved surfaces, and `transformation` is applied to the
mesh (it places the shape in the frame of its mother, as in [`place_daughter!`](@ref)).
"""
function shape_to_mesh(shape::AbstractShape{T}; facets::Integer=64,
                       transformation::Transformation3D{T}=one(Transformation3D{T}),
                       crease_angle::Real=π/6) where T
    m = GeometryBasics.mesh(shape; facets)
    pts, tris = _welded_triangles(Point3{Float64}.(coordinates(m)), faces(m))
    isempty(tris) && return GeometryBasics.Mesh(GeometryBasics.Point3{T}[], TriangleFace{Int}[];
                                                normal=GeometryBasics.Vec3f[])

    # normals of the shape at the vertices (zero for the vertices not on the surface)
    vnormals = [Vector3{Float64}(normal(shape, Point3{T}(p))) for p in pts]
    _orient!(tris, pts, vnormals)

    # unit normals of the triangles and the triangles around each vertex
    fnormals = [_triangle_normal(pts, t) for t in tris]
    vfaces = [Int[] for _ in pts]
    for (i, t) in enumerate(tris), v in t
        push!(vfaces[v], i)
    end

    # a vertex for each smooth group of triangles around a vertex
    cos_crease = cos(crease_angle)
    cos_exact = cos(crease_angle / 2)
    out_points = GeometryBasics.Point3{T}[]
    out_normals = GeometryBasics.Vec3f[]
    groups = [Tuple{Vector3{Float64},Int}[] for _ in pts]
    function vertex_index(v, f)
        nf = fnormals[f]
        n = zero(Vector3{Float64})
        for g in vfaces[v]
            fnormals[g] ⋅ nf >= cos_crease && (n += _corner_angle(pts, tris[g], v) * fnormals[g])
        end
        n = iszero(n) ? nf : normalize(n)
        for (gn, idx) in groups[v]
            gn ⋅ n > 0.999 && return idx
        end
        vn = vnormals[v]
        shading = !iszero(vn) && vn ⋅ n > cos_exact ? normalize(vn) : n
        push!(out_points, _gbpoint(Point3{T}(pts[v]) * transformation))
        push!(out_normals, GeometryBasics.Vec3f(Vector3{T}(shading) * transformation))
        push!(groups[v], (n, length(out_points)))
        return length(out_points)
    end
    out_faces = [TriangleFace{Int}(vertex_index(t[1], i), vertex_index(t[2], i), vertex_index(t[3], i))
                 for (i, t) in enumerate(tris)]
    return GeometryBasics.Mesh(out_points, out_faces; normal=out_normals)
end

# Merge the vertices closer than a small fraction of the size of the mesh, split the faces into
# triangles (fans) and drop the degenerate ones. Faces that appear twice are internal (e.g. the end
# caps of two sections of a Polycone that touch), and are dropped.
function _welded_triangles(points::AbstractVector{Point3{Float64}}, faces)
    lo = reduce((a, b) -> min.(a, b), points; init=Point3{Float64}(Inf, Inf, Inf))
    hi = reduce((a, b) -> max.(a, b), points; init=Point3{Float64}(-Inf, -Inf, -Inf))
    q = max(1e-7 * maximum(hi - lo), 1e-12)
    index = Dict{NTuple{3,Int64},Int}()
    pts = Point3{Float64}[]
    remap = map(points) do p
        get!(index, Tuple(round.(Int64, p / q))) do
            push!(pts, p)
            length(pts)
        end
    end
    welded = [[remap[Int(i)] for i in f] for f in faces]
    count = Dict{Vector{Int},Int}()
    for f in welded
        count[sort(f)] = get(count, sort(f), 0) + 1
    end
    tris = NTuple{3,Int}[]
    for f in welded
        count[sort(f)] == 1 || continue
        for k in 2:length(f)-1
            t = (f[1], f[k], f[k+1])
            allunique(t) || continue
            iszero(_triangle_cross(pts, t)) && continue
            push!(tris, t)
        end
    end
    return pts, tris
end

_triangle_cross(pts, t) = (pts[t[2]] - pts[t[1]]) × (pts[t[3]] - pts[t[1]])
_triangle_normal(pts, t) = normalize(_triangle_cross(pts, t))

function _corner_angle(pts, t, v)
    i = findfirst(==(v), t)
    a, b = pts[t[mod1(i + 1, 3)]] - pts[v], pts[t[mod1(i + 2, 3)]] - pts[v]
    return atan(norm(a × b), a ⋅ b)
end

# Orient the triangles outwards. A triangle is oriented by the normals of the shape at its vertices
# when they are on the surface. The others take the orientation of their neighbours, and the
# components with no vertex on the surface are oriented with the sign of their volume.
function _orient!(tris, pts, vnormals)
    flip(t) = (t[1], t[3], t[2])
    edges = Dict{Tuple{Int,Int},Vector{Int}}()
    for (i, t) in enumerate(tris), k in 1:3
        a, b = t[k], t[mod1(k + 1, 3)]
        push!(get!(edges, minmax(a, b), Int[]), i)
    end
    done = falses(length(tris))
    queue = Int[]
    for (i, t) in enumerate(tris)
        s = (vnormals[t[1]] + vnormals[t[2]] + vnormals[t[3]]) ⋅ _triangle_normal(pts, t)
        abs(s) > 0.1 || continue
        s < 0 && (tris[i] = flip(t))
        done[i] = true
        push!(queue, i)
    end
    # orient a neighbour so that it runs along the shared edge in the opposite direction
    function propagate!(queue)
        while !isempty(queue)
            i = popfirst!(queue)
            t = tris[i]
            for k in 1:3
                a, b = t[k], t[mod1(k + 1, 3)]
                for j in edges[minmax(a, b)]
                    done[j] && continue
                    u = tris[j]
                    any(l -> u[l] == a && u[mod1(l + 1, 3)] == b, 1:3) && (tris[j] = flip(u))
                    done[j] = true
                    push!(queue, j)
                end
            end
        end
    end
    propagate!(queue)
    for i in eachindex(tris)
        done[i] && continue
        done[i] = true
        propagate!([i])
        comp = _component(tris, edges, i)       # all the triangles connected to i
        vol = sum(t -> pts[t[1]] ⋅ (pts[t[2]] × pts[t[3]]), (tris[j] for j in comp))
        if vol < 0
            for j in comp
                tris[j] = flip(tris[j])
            end
        end
    end
    return tris
end

function _component(tris, edges, start)
    seen = Set(start)
    stack = [start]
    while !isempty(stack)
        t = tris[pop!(stack)]
        for k in 1:3, j in edges[minmax(t[k], t[mod1(k + 1, 3)])]
            j in seen || (push!(seen, j); push!(stack, j))
        end
    end
    return seen
end
