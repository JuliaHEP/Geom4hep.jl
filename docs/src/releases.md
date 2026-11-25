# Release Notes

## Unreleased

## 0.1.1 - 25-11-2026
- New features:
    - `shape_to_mesh`: mesh of a shape with normals at the vertices, used by `draw` for smooth
      shading of curved surfaces and sharp edges
      ([#24](https://github.com/JuliaHEP/Geom4hep.jl/pull/24)).
    - `csg_mesh`: mesh of the boolean shapes (and `MultiUnion`) computed with boolean operations on
      the meshes of their operands (BSP trees, as csg.js and ROOT). The result is a closed mesh,
      with the coplanar fragments merged. Used by `mesh` and `draw` for the boolean shapes and
      `MultiUnion`, which were drawn before as the overlay of their operands
      ([#25](https://github.com/JuliaHEP/Geom4hep.jl/pull/25)).
    - RoundCube example: the cube with rounded edges and corners of the Geant4.jl UserLib example,
      built from boolean shapes, with the translation of the Geant4 transformations
      ([#21](https://github.com/JuliaHEP/Geom4hep.jl/pull/21)). The Shapes and RoundCube examples
      draw the boolean shapes ([#25](https://github.com/JuliaHEP/Geom4hep.jl/pull/25)).
- Changes:
    - `draw` and `draw!` pass the extra keyword arguments to the Makie plotting functions, e.g.
      `draw(shape; color=:blue)` or `draw(shape; wireframe=true, linewidth=2)`
      ([#26](https://github.com/JuliaHEP/Geom4hep.jl/pull/26)).
- Fixes:
    - Boolean shapes: `inside` and `normal` at the points on the surface of both operands (seams of
      unions, edges where the cut of a subtraction meets the surface), and the quick bounding box
      checks for the points on the surface of an operand
      ([#21](https://github.com/JuliaHEP/Geom4hep.jl/pull/21)).
- Maintenance:
    - `[compat]` entries for the docs environment (CompatHelper:
      [#18](https://github.com/JuliaHEP/Geom4hep.jl/pull/18),
      [#19](https://github.com/JuliaHEP/Geom4hep.jl/pull/19),
      [#20](https://github.com/JuliaHEP/Geom4hep.jl/pull/20),
      [#23](https://github.com/JuliaHEP/Geom4hep.jl/pull/23)).

## 0.1.0 - 19-11-2026
- New features:
    - Shapes: `Sphere` (with ϕ and θ segments, also GDML `orb`)
      ([#9](https://github.com/JuliaHEP/Geom4hep.jl/pull/9)), `Polyhedra`
      ([#10](https://github.com/JuliaHEP/Geom4hep.jl/pull/10)) and `Torus`
      ([#12](https://github.com/JuliaHEP/Geom4hep.jl/pull/12)), and the mesh of all the shapes for
      drawing.
    - Flat geometry model (`close_geometry`) and navigator (`locate`, `compute_step`) shared by CPUs
      and GPUs (KernelAbstractions), with boolean shapes as flat CSG programs and BVHs for the
      daughters of volumes. `xray` images on CPUs and GPUs
      ([#8](https://github.com/JuliaHEP/Geom4hep.jl/pull/8)).
    - Single precision (`Float32`) support, needed by Apple GPUs
      ([#7](https://github.com/JuliaHEP/Geom4hep.jl/pull/7),
      [#8](https://github.com/JuliaHEP/Geom4hep.jl/pull/8)).
    - Shape benchmarks with the methodology of VecGeom (`benchmark_shape`, `vecgeom_suite`) and
      comparison with VecGeom ([#14](https://github.com/JuliaHEP/Geom4hep.jl/pull/14)).
    - Documentation with Documenter and Literate examples
      ([#15](https://github.com/JuliaHEP/Geom4hep.jl/pull/15)).
- Changes:
    - `snake_case` API (`distance_to_in`, `read_gdml`, `world_volume`, ...) and own point and vector
      types (`Point3`, `Vector3`) ([#7](https://github.com/JuliaHEP/Geom4hep.jl/pull/7)).
    - Drawing functions moved to a Makie package extension
      ([#7](https://github.com/JuliaHEP/Geom4hep.jl/pull/7)).
    - Repository moved to the JuliaHEP organization
      ([#16](https://github.com/JuliaHEP/Geom4hep.jl/pull/16)), and prepared for registration in
      the General registry: license, TagBot and CompatHelper
      ([#17](https://github.com/JuliaHEP/Geom4hep.jl/pull/17)).
- Fixes:
    - Many fixes of the shape functions at surfaces and for grazing rays
      ([#7](https://github.com/JuliaHEP/Geom4hep.jl/pull/7),
      [#8](https://github.com/JuliaHEP/Geom4hep.jl/pull/8),
      [#9](https://github.com/JuliaHEP/Geom4hep.jl/pull/9)).
    - Navigator: the X-ray kernel of the full CMS geometry compiles again on Metal
      ([#11](https://github.com/JuliaHEP/Geom4hep.jl/pull/11)).
    - `Trap` with a face that has an edge of zero length (CMS TEC module frames)
      ([#13](https://github.com/JuliaHEP/Geom4hep.jl/pull/13)).
