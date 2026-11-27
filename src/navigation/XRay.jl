#---X-ray image of a geometry (benchmark of the navigation, as the XRay benchmark of VecGeom)-----------
#   Rays parallel to one of the axes cross the world volume; each pixel accumulates the length of the
#   traversed volumes weighted by a quantity of their material: the mass thickness (density × length),
#   or the material budget in radiation lengths (length/X0) or interaction lengths (length/λ). The same per-pixel function runs in a
#   KernelAbstractions kernel on GPUs, or in a loop over Julia threads on the CPU.
#
#   Note: the CPU backend of KernelAbstractions ≥ 0.10 (POCL) inlines all functions, which makes the
#   kernels explode in size for geometries with nested boolean shapes; Julia threads are used instead.

@inline function _xray_pixel(model::GeometryModel{T}, weights, origin, du, dv, dir, maxsteps, i, j) where T
    point = origin + T(i - 1) * du + T(j - 1) * dv
    state = locate(model, point)
    sum = zero(T)
    n = Int32(0)
    while in_world(state) && n < maxsteps
        w = @inbounds weights[material_index(model, state)]
        step, state = compute_step(model, state, point, dir, T(Inf))
        point = point + step * dir
        sum += step * w
        n += Int32(1)
    end
    return sum, n
end

# Weight of the length in each material of the model for the quantity of the X-ray image
function _xray_weights(model::GeometryModel{T}, quantity::Symbol) where T
    mats = Array(model.materials)
    quantity == :mass   && return T[m.density for m in mats]
    quantity == :radlen && return T[1 / m.radlen for m in mats]
    quantity == :nuclen && return T[1 / m.nuclen for m in mats]
    throw(ArgumentError("unknown X-ray quantity :$quantity (use :mass, :radlen or :nuclen)"))
end

# The image is processed in tiles of _TILE_I × _TILE_J pixels (adjacent pixels, similar rays), visited
# in a scrambled order (see `xray`)
const _TILE_I = 64
const _TILE_J = 16

@kernel function _xray_kernel!(image, nsteps, model, weights, origin, du, dv, dir, maxsteps, t0, ntiles, stride)
    g = @index(Global, Linear) - 1
    k, w = divrem(g, _TILE_I * _TILE_J)                 # tile in the batch, pixel in the tile
    c, r = divrem(w, _TILE_I)
    tile = mod((k + t0) * Int64(stride), Int64(ntiles))  # tile in the image
    tj, ti = divrem(Int(tile), cld(size(image, 1), _TILE_I))
    i, j = ti * _TILE_I + r + 1, tj * _TILE_J + c + 1
    if i <= size(image, 1) && j <= size(image, 2)
        sum, n = _xray_pixel(model, weights, origin, du, dv, dir, maxsteps, i, j)
        @inbounds image[i, j] = sum
        @inbounds nsteps[i, j] = n
    end
end

"""
    xray(model::GeometryModel, axis, npixels; quantity=:mass, backend=CPU(), maxsteps=10_000, launch_time=0.5)

X-ray image of the world volume of `model` with rays parallel to `axis` (1, 2 or 3 for x, y, z):
the `quantity` integrated along each ray,
- `:mass`: the mass thickness (density × length, in g/cm³ × mm);
- `:radlen`: the material budget in radiation lengths (length/X0);
- `:nuclen`: the material budget in nuclear interaction lengths (length/λ).
The materials without composition (only a density) have unknown X0 and λ, which make the pixels
that cross them `NaN` (the ones with zero density, e.g. a vacuum, do not contribute).
`npixels` is the number of pixels
along the first transverse axis. Returns `(u, v, image, nsteps)` with the pixel centres along the two
transverse axes, the image and the number of navigation steps per pixel.

`backend` is a KernelAbstractions backend: `CPU()` (Julia threads), or a GPU backend such as
`MetalBackend()` or `CUDABackend()` (the model is copied to the device). On GPUs the image is
computed in several kernel launches, each one taking about `launch_time` seconds: a long-running
kernel on a GPU that also drives a display is aborted by the operating system (e.g. macOS error
"Impacting Interactivity").
"""
function xray(model::GeometryModel{T}, axis::Integer, npixels::Integer; quantity::Symbol=:mass, backend=CPU(),
              maxsteps::Integer=10_000, launch_time::Real=0.5) where T
    weights = _xray_weights(model, quantity)
    lower, upper = shape_call(extent, model, model.volumes[1].shape)
    iu, iv = axis == 1 ? (2, 3) : axis == 2 ? (1, 3) : (1, 2)
    nu = npixels
    nv = max(1, round(Int, nu * (upper[iv] - lower[iv]) / (upper[iu] - lower[iu])))
    δu, δv = (upper[iu] - lower[iu]) / nu, (upper[iv] - lower[iv]) / nv
    unit(k) = Vector3{T}(k == 1, k == 2, k == 3)
    origin = lower + push_tolerance(T) * unit(axis) + (δu / 2) * unit(iu) + (δv / 2) * unit(iv)
    du, dv, dir = δu * unit(iu), δv * unit(iv), unit(axis)
    if backend isa CPU
        image = zeros(T, nu, nv)
        nsteps = zeros(Int32, nu, nv)
        Threads.@threads for j in 1:nv
            for i in 1:nu
                image[i, j], nsteps[i, j] = _xray_pixel(model, weights, origin, du, dv, dir, Int32(maxsteps), i, j)
            end
        end
    else
        dimage = KernelAbstractions.zeros(backend, T, nu, nv)
        dnsteps = KernelAbstractions.zeros(backend, Int32, nu, nv)
        dmodel = adapt(backend, model)
        dweights = adapt(backend, weights)
        kernel! = _xray_kernel!(backend)
        # The tiles are processed in batches whose size adapts to take ~launch_time each. The cost of
        # a ray varies a lot across the image (by orders of magnitude between the empty border and the
        # centre of a detector), so the tiles are visited in a scrambled order (k ↦ k⋅stride mod
        # ntiles, stride coprime with ntiles): every batch is a sample of the whole image, with a
        # predictable cost.
        ntiles = cld(nu, _TILE_I) * cld(nv, _TILE_J)
        stride = max(1, round(Int, 0.618 * ntiles))
        while gcd(stride, ntiles) != 1
            stride += 1
        end
        nbatch, t0 = 1, 0
        while t0 < ntiles
            n = min(nbatch, ntiles - t0)
            t = @elapsed begin
                kernel!(dimage, dnsteps, dmodel, dweights, origin, du, dv, dir, Int32(maxsteps), t0, ntiles, stride;
                        ndrange=n * _TILE_I * _TILE_J)
                KernelAbstractions.synchronize(backend)
            end
            t0 += n
            nbatch = clamp(round(Int, n * launch_time / max(t, 1e-3)), 1, 2n)   # grow at most ×2
        end
        image, nsteps = Array(dimage), Array(dnsteps)
    end
    u = range(lower[iu] + δu / 2, step=δu, length=nu)
    v = range(lower[iv] + δv / 2, step=δv, length=nv)
    return u, v, image, nsteps
end
