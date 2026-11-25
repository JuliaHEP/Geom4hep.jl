#---Boolean operations on meshes with BSP trees-----------------------------------------------------
#   A port of the csg.js algorithm (Evan Wallace, MIT license), as used by ROOT (RootCsg) to draw
#   composite shapes. A solid is a set of convex polygons, oriented outwards, that enclose it. A BSP
#   tree of the polygons of a solid A is used to clip away the polygons of a solid B that are inside
#   (or outside) A. All operations are combinations of clipping, inverting (swapping inside and
#   outside) and building trees:
#     A ∪ B : A - polygons of A inside B, plus the polygons of B outside A
#     A - B : ~(~A ∪ B)
#     A ∩ B : ~(~A ∪ ~B)
#   The coplanar polygons are classified with a tolerance `eps` (the distance to the plane), and
#   put in front or behind according to the orientation of their normal. Polygons split by a plane
#   keep the plane of the original polygon.

struct CSGPlane
    n::Vector3{Float64}
    w::Float64
end
Base.:-(p::CSGPlane) = CSGPlane(-p.n, -p.w)

struct CSGPolygon
    vertices::Vector{Point3{Float64}}
    plane::CSGPlane
end
function CSGPolygon(vertices::Vector{Point3{Float64}})
    a, b, c = vertices[1], vertices[2], vertices[3]
    n = normalize((b - a) × (c - a))
    return CSGPolygon(vertices, CSGPlane(n, n ⋅ a))
end
flip(p::CSGPolygon) = CSGPolygon(reverse(p.vertices), -p.plane)

const COPLANAR, FRONT, BACK, SPANNING = 0, 1, 2, 3

# Split `poly` by `plane`, and put the pieces in the lists of polygons in front or behind it. The
# coplanar polygons go into `cofront` or `coback` according to their orientation.
function split_polygon!(plane::CSGPlane, poly::CSGPolygon, cofront, coback, front, back, eps)
    vs = poly.vertices
    nv = length(vs)
    types = Vector{Int}(undef, nv)
    ptype = COPLANAR
    for (i, v) in enumerate(vs)
        t = plane.n ⋅ v - plane.w
        types[i] = t < -eps ? BACK : t > eps ? FRONT : COPLANAR
        ptype |= types[i]
    end
    if ptype == COPLANAR
        push!(plane.n ⋅ poly.plane.n > 0 ? cofront : coback, poly)
    elseif ptype == FRONT
        push!(front, poly)
    elseif ptype == BACK
        push!(back, poly)
    else
        f, b = Point3{Float64}[], Point3{Float64}[]
        for i in 1:nv
            j = mod1(i + 1, nv)
            ti, tj, vi, vj = types[i], types[j], vs[i], vs[j]
            ti != BACK && push!(f, vi)
            ti != FRONT && push!(b, vi)
            if ti | tj == SPANNING
                t = (plane.w - plane.n ⋅ vi) / (plane.n ⋅ (vj - vi))
                v = vi + (vj - vi) * t
                push!(f, v)
                push!(b, v)
            end
        end
        length(f) >= 3 && push!(front, CSGPolygon(f, poly.plane))
        length(b) >= 3 && push!(back, CSGPolygon(b, poly.plane))
    end
    return nothing
end

mutable struct BSPNode
    plane::Union{Nothing,CSGPlane}
    front::Union{Nothing,BSPNode}
    back::Union{Nothing,BSPNode}
    polygons::Vector{CSGPolygon}
end
BSPNode() = BSPNode(nothing, nothing, nothing, CSGPolygon[])
function BSPNode(polygons::Vector{CSGPolygon}, eps)
    node = BSPNode()
    build!(node, polygons, eps)
    return node
end

# All the nodes of the tree (iterative, the trees of convex solids are as deep as their polygons)
function nodes(root::BSPNode)
    result = BSPNode[]
    stack = [root]
    while !isempty(stack)
        node = pop!(stack)
        push!(result, node)
        isnothing(node.front) || push!(stack, node.front)
        isnothing(node.back) || push!(stack, node.back)
    end
    return result
end

all_polygons(root::BSPNode) = reduce(vcat, (node.polygons for node in nodes(root)); init=CSGPolygon[])

# Swap the inside and the outside of the solid
function invert!(root::BSPNode)
    for node in nodes(root)
        node.polygons = flip.(node.polygons)
        isnothing(node.plane) || (node.plane = -node.plane)
        node.front, node.back = node.back, node.front
    end
    return root
end

# The parts of `polygons` that are outside the solid of the tree
function clip_polygons(root::BSPNode, polygons::Vector{CSGPolygon}, eps)
    result = CSGPolygon[]
    stack = [(root, polygons)]
    while !isempty(stack)
        node, polys = pop!(stack)
        if isnothing(node.plane)
            append!(result, polys)
            continue
        end
        front, back = CSGPolygon[], CSGPolygon[]
        for p in polys
            split_polygon!(node.plane, p, front, back, front, back, eps)
        end
        if isnothing(node.front)
            append!(result, front)
        else
            push!(stack, (node.front, front))
        end
        isnothing(node.back) || push!(stack, (node.back, back))    # behind a leaf: inside, dropped
    end
    return result
end

# Remove the polygons of the tree `root` that are inside the solid of `other`
function clip_to!(root::BSPNode, other::BSPNode, eps)
    for node in nodes(root)
        node.polygons = clip_polygons(other, node.polygons, eps)
    end
    return root
end

function build!(root::BSPNode, polygons::Vector{CSGPolygon}, eps)
    stack = [(root, polygons)]
    while !isempty(stack)
        node, polys = pop!(stack)
        isempty(polys) && continue
        isnothing(node.plane) && (node.plane = polys[1].plane)
        front, back = CSGPolygon[], CSGPolygon[]
        for p in polys
            split_polygon!(node.plane, p, node.polygons, node.polygons, front, back, eps)
        end
        if !isempty(front)
            isnothing(node.front) && (node.front = BSPNode())
            push!(stack, (node.front, front))
        end
        if !isempty(back)
            isnothing(node.back) && (node.back = BSPNode())
            push!(stack, (node.back, back))
        end
    end
    return root
end

function csg_union(a::Vector{CSGPolygon}, b::Vector{CSGPolygon}, eps)
    A, B = BSPNode(a, eps), BSPNode(b, eps)
    clip_to!(A, B, eps)
    clip_to!(B, A, eps)
    invert!(B); clip_to!(B, A, eps); invert!(B)
    build!(A, all_polygons(B), eps)
    return all_polygons(A)
end

function csg_subtraction(a::Vector{CSGPolygon}, b::Vector{CSGPolygon}, eps)
    A, B = BSPNode(a, eps), BSPNode(b, eps)
    invert!(A)
    clip_to!(A, B, eps)
    clip_to!(B, A, eps)
    invert!(B); clip_to!(B, A, eps); invert!(B)
    build!(A, all_polygons(B), eps)
    invert!(A)
    return all_polygons(A)
end

function csg_intersection(a::Vector{CSGPolygon}, b::Vector{CSGPolygon}, eps)
    A, B = BSPNode(a, eps), BSPNode(b, eps)
    invert!(A)
    clip_to!(B, A, eps)
    invert!(B)
    clip_to!(A, B, eps)
    clip_to!(B, A, eps)
    build!(A, all_polygons(B), eps)
    invert!(A)
    return all_polygons(A)
end

#---Polygons of shapes-------------------------------------------------------------------------------
# The triangles of the mesh of a primitive shape, welded and oriented outwards
function csg_polygons(shape::AbstractShape{T}, t::Transformation3D{T}, facets, eps) where T
    m = GeometryBasics.mesh(shape; facets)
    pts, tris = _welded_triangles(Point3{Float64}.(coordinates(m)), faces(m))
    vnormals = [Vector3{Float64}(normal(shape, Point3{T}(p))) for p in pts]
    _orient!(tris, pts, vnormals)
    tpts = [Point3{Float64}(Point3{T}(p) * t) for p in pts]
    return [CSGPolygon([tpts[i] for i in tri]) for tri in tris]
end

function csg_polygons(shape::BooleanShape{T,OP}, t::Transformation3D{T}, facets, eps) where {T,OP}
    a = csg_polygons(shape.left, t, facets, eps)
    b = csg_polygons(shape.right, shape.transformation * t, facets, eps)
    OP === :union && return csg_union(a, b, eps)
    OP === :subtraction && return csg_subtraction(a, b, eps)
    return csg_intersection(a, b, eps)
end

function csg_polygons(mu::MultiUnion{T}, t::Transformation3D{T}, facets, eps) where T
    isempty(mu.shapes) && return CSGPolygon[]
    polys = csg_polygons(mu.shapes[1], mu.transformations[1] * t, facets, eps)
    for i in 2:length(mu.shapes)
        polys = csg_union(polys, csg_polygons(mu.shapes[i], mu.transformations[i] * t, facets, eps), eps)
    end
    return polys
end

"""
    csg_mesh(shape; facets=64, eps=1e-7)

Surface mesh of a boolean shape (or a [`MultiUnion`](@ref)) computed with boolean operations on the
meshes of its operands (BSP trees, as csg.js). `eps` is the tolerance to classify vertices on the
planes of the polygons, relative to the size of the shape.

The polygons of the result are cleaned up:
- the vertices are welded, and the T-junctions (vertices of a polygon in the middle of an edge of
  its neighbour, left by the splits) are removed by inserting these vertices in the edges;
- the connected coplanar fragments are merged into faces (with holes), and the vertices in the
  middle of straight edges are removed, and the faces are triangulated again.
The mesh is closed: every edge is shared by two triangles, in opposite directions.
"""
function csg_mesh(shape::AbstractShape{T}; facets::Integer=64, eps::Real=1e-7) where T
    lo, hi = extent(shape)
    abseps = eps * max(maximum(abs.(hi - lo)), 1e-12)
    polys = csg_polygons(shape, one(Transformation3D{T}), facets, abseps)
    pts, tris = _clean_polygons(polys, abseps)
    return GeometryBasics.Mesh(_gbpoint.(Point3{T}.(pts)), tris)
end

#---Cleanup of the polygons of the result---------------------------------------------------------
# The vertices and the triangles of the closed mesh made of the polygons `polys` (see csg_mesh)
function _clean_polygons(polys::Vector{CSGPolygon}, tol)
    pts, loops, planes = _weld_polygons(polys, tol)
    _drop_degenerate!(loops, planes, pts, tol)
    _insert_tjunctions!(loops, pts, tol)
    faces, fplanes = _merge_coplanar(loops, planes, pts, tol)
    _remove_collinear!(faces, pts, tol)
    loops = reduce(vcat, (_face_pieces(f, p, pts, tol) for (f, p) in zip(faces, fplanes)); init=Vector{Int}[])
    return pts, _triangulate_loops!(loops, pts, tol)
end

# Triangles of convex loops (the pieces of the faces) that make a closed surface. The slivers of
# zero width (e.g. from earcut) are dropped, and their vertices are inserted in the edges of the
# neighbouring loops (with the vertices skipped by earcut), before the loops are triangulated.
function _triangulate_loops!(loops, pts, tol)
    _drop_degenerate!(loops, nothing, pts, tol)
    _insert_tjunctions!(loops, pts, tol)
    tris = TriangleFace{Int}[]
    for loop in loops
        _triangulate_convex!(tris, pts, loop, tol)
    end
    return tris
end

# Merge the vertices closer than `tol` (a hash grid of cells of size `tol`, looking in the 27
# neighbouring cells). The polygons become loops of indices, without repeated vertices.
function _weld_polygons(polys::Vector{CSGPolygon}, tol)
    pts = Point3{Float64}[]
    grid = Dict{NTuple{3,Int},Vector{Int}}()
    cell(p) = (floor(Int, p[1] / tol), floor(Int, p[2] / tol), floor(Int, p[3] / tol))
    function index(p)
        c = cell(p)
        for dx in -1:1, dy in -1:1, dz in -1:1
            for i in get(grid, (c[1] + dx, c[2] + dy, c[3] + dz), ())
                norm(pts[i] - p) <= tol && return i
            end
        end
        push!(pts, p)
        push!(get!(grid, c, Int[]), length(pts))
        return length(pts)
    end
    loops = Vector{Int}[]
    planes = CSGPlane[]
    for poly in polys
        loop = Int[]
        for v in poly.vertices
            i = index(v)
            (isempty(loop) || loop[end] != i) && push!(loop, i)
        end
        length(loop) > 1 && loop[end] == loop[1] && pop!(loop)
        length(loop) >= 3 && (push!(loops, loop); push!(planes, poly.plane))
    end
    return pts, loops, planes
end

# Remove the loops with all their vertices within `tol` of a line (slivers of zero width). Their
# vertices are then inserted in the edges of the neighbours by _insert_tjunctions!.
function _drop_degenerate!(loops, planes, pts, tol)
    keep = map(loops) do loop
        a = norm(sum(e -> pts[e[1]] × pts[e[2]], _loop_edges(loop); init=zero(Vector3{Float64})))
        longest = maximum(e -> norm(pts[e[2]] - pts[e[1]]), _loop_edges(loop))
        a > tol * longest                    # twice the area > tolerance * length
    end
    isnothing(planes) || keepat!(planes, keep)
    keepat!(loops, keep)
    return loops
end

_loop_edges(loop) = ((loop[k], loop[mod1(k + 1, length(loop))]) for k in eachindex(loop))

# Signed area of a loop, positive when it turns counterclockwise around the normal `n`
_loop_area(loop, pts, n) = sum(e -> pts[e[1]] × pts[e[2]], _loop_edges(loop); init=zero(Vector3{Float64})) ⋅ n / 2

# Merge the connected polygons on the same plane into faces: a list of loops, the outer boundary
# first and then the holes. A group of polygons is left unmerged when its boundary is ambiguous
# (two loops touching at a vertex) or it has not exactly one outer boundary.
function _merge_coplanar(loops::Vector{Vector{Int}}, planes::Vector{CSGPlane}, pts, tol)
    n = length(loops)
    parent = collect(1:n)
    function root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    owner = Dict{Tuple{Int,Int},Int}()
    for (l, loop) in enumerate(loops), e in _loop_edges(loop)
        owner[e] = l
    end
    for ((a, b), l) in owner
        m = get(owner, (b, a), 0)
        m > 0 || continue
        p, q = planes[l], planes[m]
        p.n ⋅ q.n > 1 - 1e-9 && abs(p.w - q.w) <= tol && (parent[root(l)] = root(m))
    end
    groups = Dict{Int,Vector{Int}}()
    for l in 1:n
        push!(get!(groups, root(l), Int[]), l)
    end
    faces, fplanes = Vector{Vector{Int}}[], CSGPlane[]
    for l in sort!([members[1] for members in values(groups)])     # deterministic order
        members = groups[root(l)]
        plane = planes[members[1]]
        boundary = length(members) > 1 ? _group_boundary(members, loops, owner, root, plane, pts) : nothing
        if isnothing(boundary)
            for m in members
                push!(faces, [loops[m]]); push!(fplanes, planes[m])
            end
        else
            push!(faces, boundary); push!(fplanes, plane)
        end
    end
    return faces, fplanes
end

# The boundary loops of a group of coplanar polygons (the edges not shared within the group)
function _group_boundary(members, loops, owner, root, plane, pts)
    g = root(members[1])
    next = Dict{Int,Int}()
    for m in members, (a, b) in _loop_edges(loops[m])
        t = get(owner, (b, a), 0)
        t > 0 && root(t) == g && continue
        haskey(next, a) && return nothing                 # two boundary edges leave a: ambiguous
        next[a] = b
    end
    result = Vector{Int}[]
    while !isempty(next)
        start = first(keys(next))
        loop = [start]
        v = pop!(next, start)
        while v != start
            haskey(next, v) || return nothing
            push!(loop, v)
            v = pop!(next, v)
        end
        push!(result, loop)
    end
    areas = [_loop_area(loop, pts, plane.n) for loop in result]
    count(>(0), areas) == 1 || return nothing
    return result[sortperm(areas, rev=true)]                # the outer boundary first
end

# Remove the vertices in the middle of straight edges. A vertex is removed only when it is in the
# middle of a straight edge in all the loops it belongs to, so that the neighbouring faces agree.
function _remove_collinear!(faces, pts, tol)
    removable = Dict{Int,Bool}()
    for face in faces, loop in face
        n = length(loop)
        for k in 1:n
            a, v, b = pts[loop[mod1(k - 1, n)]], pts[loop[k]], pts[loop[mod1(k + 1, n)]]
            straight = norm((v - a) × (b - v)) <= tol * norm(b - a) && (v - a) ⋅ (b - v) > 0
            removable[loop[k]] = get(removable, loop[k], true) && straight
        end
    end
    for face in faces, loop in face
        filter!(v -> !removable[v], loop)
    end
    for face in faces
        filter!(loop -> length(loop) >= 3, face)
    end
    filter!(!isempty, faces)
    return faces
end

# Convex pieces of a face (loops of vertices). A face with a single convex loop is kept, the others
# are triangulated with earcut, in the projection on the coordinate plane closest to the face.
function _face_pieces(face, plane::CSGPlane, pts, tol)
    isempty(face) && return Vector{Int}[]
    if length(face) == 1
        loop = face[1]
        n = length(loop)
        convex = all(1:n) do k
            a, v, b = pts[loop[mod1(k - 1, n)]], pts[loop[k]], pts[loop[mod1(k + 1, n)]]
            ((v - a) × (b - v)) ⋅ plane.n >= -tol * norm(b - a)
        end
        convex && return [loop]
    end
    ax = argmax(abs.(plane.n))
    i, j = ax == 1 ? (2, 3) : ax == 2 ? (3, 1) : (1, 2)
    flat = reduce(vcat, face)
    poly = [[GeometryBasics.Point{2,Float64}(pts[v][i], pts[v][j]) for v in loop] for loop in face]
    pieces = Vector{Int}[]
    for t in GeometryBasics.earcut_triangulate(poly)
        a, b, c = (flat[convert(Int, k)] for k in t)
        ((pts[b] - pts[a]) × (pts[c] - pts[a])) ⋅ plane.n < 0 && ((b, c) = (c, b))
        push!(pieces, [a, b, c])
    end
    return pieces
end

# Insert in the edges without a twin (the same edge, in the opposite direction, in another polygon)
# the vertices that lie on them, ordered along the edge.
function _insert_tjunctions!(loops, pts::Vector{Point3{Float64}}, tol)
    edges = Set{Tuple{Int,Int}}()
    for loop in loops, k in eachindex(loop)
        push!(edges, (loop[k], loop[mod1(k + 1, length(loop))]))
    end
    used = sort!(unique!(reduce(vcat, loops; init=Int[])), by=i -> pts[i][1])
    xs = [pts[i][1] for i in used]
    for (l, loop) in enumerate(loops)
        any(k -> (loop[mod1(k + 1, length(loop))], loop[k]) ∉ edges, eachindex(loop)) || continue
        newloop = Int[]
        for k in eachindex(loop)
            a, b = loop[k], loop[mod1(k + 1, length(loop))]
            push!(newloop, a)
            (b, a) in edges && continue
            pa, pb = pts[a], pts[b]
            d = pb - pa
            len2 = d ⋅ d
            xlo, xhi = minmax(pa[1], pb[1])
            onedge = Tuple{Float64,Int}[]
            for j in searchsortedfirst(xs, xlo - tol):searchsortedlast(xs, xhi + tol)
                i = used[j]
                (i == a || i == b) && continue
                t = (pts[i] - pa) ⋅ d / len2
                0 < t < 1 || continue
                norm(pts[i] - (pa + d * t)) <= tol && push!(onedge, (t, i))
            end
            sort!(onedge)
            append!(newloop, last.(onedge))
        end
        loops[l] = newloop
    end
    return loops
end

# Triangles of a convex polygon that may have collinear vertices (the inserted T-junctions). Ears
# are cut at the corners, choosing the best shaped triangle (largest smallest angle) that leaves a
# polygon with a non-zero area. When there is no such ear, the polygon is fanned from its centroid.
function _triangulate_convex!(tris, pts, loop, tol)
    loop = copy(loop)
    area(l) = norm(sum(k -> pts[l[k]] × pts[l[mod1(k + 1, length(l))]], eachindex(l); init=zero(Vector3{Float64}))) / 2
    function minangle(a, b, c)
        u, v, w = pts[b] - pts[a], pts[c] - pts[b], pts[a] - pts[c]
        s = norm(u × v)
        return min(atan(s, -(u ⋅ w)), atan(s, -(u ⋅ v)), atan(s, -(v ⋅ w)))
    end
    while length(loop) > 3
        n = length(loop)
        best, bestq = 0, -1.0
        for k in 1:n
            a, b, c = loop[mod1(k - 1, n)], loop[k], loop[mod1(k + 1, n)]
            norm((pts[b] - pts[a]) × (pts[c] - pts[b])) > tol * norm(pts[c] - pts[a]) || continue
            q = minangle(a, b, c)
            q > bestq || continue
            rest = deleteat!(copy(loop), k)
            area(rest) > tol * sqrt(area(loop)) || continue
            best, bestq = k, q
        end
        if best == 0                                    # no valid ear: fan from the centroid
            push!(pts, sum(pts[i] for i in loop) / n)
            c = length(pts)
            for k in 1:n
                push!(tris, TriangleFace{Int}(loop[k], loop[mod1(k + 1, n)], c))
            end
            return tris
        end
        push!(tris, TriangleFace{Int}(loop[mod1(best - 1, n)], loop[best], loop[mod1(best + 1, n)]))
        deleteat!(loop, best)
    end
    area(loop) > 0 && push!(tris, TriangleFace{Int}(loop[1], loop[2], loop[3]))
    return tris
end
