module Geom4hepMakieExt

using Geom4hep
using Geom4hep: AbstractShape, Volume, Transformation3D, AABB, BVH, build_bvh, extent, inside,
                kInside, kOutside, distance_to_in, distance_to_out, Point3, Vector3, _gbpoint,
                shape_to_mesh
import Geom4hep: draw, draw!, draw_distance_to_in, draw_distance_to_out
using GeometryBasics: GeometryBasics, Mesh, coordinates, faces
using LinearAlgebra: normalize
using Makie

# Grey levels from light (outer volumes) to dark (deep daughters)
const colors = [Makie.RGBf(g, g, g) for g in range(0.95f0, 0.05f0, length=16)]

#---Draw a Volume---------------------------------------------------------------
function draw!(s::LScene, vol::Volume{T}, t::Transformation3D{T}, level::Int64, wireframe::Bool, bvh::Bool, maxlevel::Int64; kwargs...) where T
    if wireframe                       # the faces of the shape, without the diagonals of triangles
        m = GeometryBasics.mesh(vol.shape)
        isone(t) || (m = Mesh(map(c -> _gbpoint(Point3{T}(c) * t), coordinates(m)), faces(m)))
    else
        m = shape_to_mesh(vol.shape; transformation=t)
    end
    icolor = clamp(level, 1, length(colors))
    if wireframe
        wireframe!(s, m; color=colors[end+1-icolor], visible = level != 1, kwargs...)
    else
        mesh!(s, m; color=colors[icolor], transparency=true, visible = level != 1, kwargs...)
    end
    if level < maxlevel
        for daughter in vol.daughters
            if level == 1 && bvh
                draw!(s, Geom4hep.AABB(daughter))
            end
            draw!(s, daughter.volume, daughter.transformation * t, level+1, wireframe, bvh, maxlevel; kwargs...)
        end
        if level == 1 && bvh && !isempty(vol.daughters)
            draw!(s, build_bvh(vol.daughters))
        end
    end
    return s
end

function draw!(s::LScene, vol::Volume{T}; wireframe::Bool=false, bvh::Bool=false, maxlevel::Int64=999, kwargs...) where T
    draw!(s, vol, one(Transformation3D{T}), 1, wireframe, bvh, maxlevel; kwargs...)
end

function draw(vol::Volume{T}; wireframe::Bool=false, bvh::Bool=false, maxlevel::Int64=999, kwargs...) where T
    fig = Figure(size=(1280, 720))
    s = LScene(fig[1,1])
    draw!(s, vol, one(Transformation3D{T}), 1, wireframe, bvh, maxlevel; kwargs...)
    return fig
end

#---Draw a Shape---------------------------------------------------------------
function draw!(s::LScene, shape::AbstractShape; wireframe::Bool=false, kwargs...)
    m = wireframe ? GeometryBasics.mesh(shape) : shape_to_mesh(shape)
    if wireframe
        wireframe!(s, m; kwargs...)
    else
        mesh!(s, m; color=colors[2], transparency=true, kwargs...)
    end
    return s
end

function draw(shape::AbstractShape; wireframe::Bool=false, kwargs...)
    fig = Figure()
    s = LScene(fig[1,1])
    draw!(s, shape; wireframe, kwargs...)
    return fig
end

#---Draw a AABB------------------------------------------------------------------
function draw!(s::LScene, aabb::AABB; color=:blue, kwargs...)
    wireframe!(s, Mesh(_gbpoint.(coordinates(aabb)), faces(aabb)); color, kwargs...)
    return s
end

function draw!(s::LScene, bvh::BVH; color=:green, kwargs...)
    if isempty(bvh.children)
        draw!(s, bvh.aabb; color, kwargs...)
    else
        for c in bvh.children
            draw!(s, c; color, kwargs...)
        end
    end
    return s
end

function draw(bvh::BVH; color=:green, kwargs...)
    fig = Figure()
    s = LScene(fig[1,1])
    draw!(s, bvh; color, kwargs...)
    return fig
end

#---Testing functions------------------------------------------------------------
function draw_distance_to_out(shape::AbstractShape{T}, N::Integer) where T<:AbstractFloat
    low, hi = extent(shape)
    dim = hi - low
    points = (low + rp * dim for rp in rand(Vector3{T}, N))
    result = Vector{Point3{T}}()
    for point in points
        if inside(shape, point) == kInside
            dir = normalize(rand(Vector3{T}) + Vector3{T}(-.5,-.5,-.5))
            push!(result, point + dir * distance_to_out(shape, point, dir))
        end
    end
    fig = Figure()
    s = LScene(fig[1, 1])
    scatter!(s, _gbpoint.(result), color=:black, markersize=1)
    scatter!(s, _gbpoint.([low, hi]), color=:blue, markersize=10)
    return fig
end

function draw_distance_to_in(shape::AbstractShape{T}, N::Integer) where T<:AbstractFloat
    low, hi = extent(shape)
    dim = hi - low
    low -= dim/10
    hi  += dim/10
    dim = hi - low
    points = (low + rp * dim for rp in rand(Vector3{T}, N))
    dirs = normalize.(rand(Vector3{T}, N) .+ Ref(Vector3{T}(-.5,-.5,-.5)))
    result = Vector{Point3{T}}()
    for (point, dir) in zip(points, dirs)
        if inside(shape, point) == kOutside
            dist = distance_to_in(shape, point, dir)
            if dist != Inf
                push!(result, point + dir * dist)
            end
        end
    end
    fig = Figure()
    s = LScene(fig[1, 1])
    scatter!(s, _gbpoint.(result), color=:black, markersize=1)
    scatter!(s, _gbpoint.([extent(shape)...]), color=:blue, markersize=10)
    return fig
end

end # module
