using LinearAlgebra
import GeometryBasics

@testset "shape_to_mesh{$T}" for T in (Float64, Float32)
    signed_volume(P, F) = sum(f -> P[f[1]] ⋅ (P[f[2]] × P[f[3]]), F) / 6
    @testset "$(shapename(s))" for s in test_shapes(T)
        m = shape_to_mesh(s)
        P, F, N = GeometryBasics.coordinates(m), GeometryBasics.faces(m), m.normal
        @test !isempty(F)
        @test length(N) == length(P)
        @test all(n -> norm(n) ≈ 1, N)
        # the normals of the vertices are on the outer side of their triangles
        @test all(f -> ((P[f[2]] - P[f[1]]) × (P[f[3]] - P[f[1]])) ⋅ (N[f[1]] + N[f[2]] + N[f[3]]) > 0, F)
        # a closed surface oriented outwards encloses the volume of the shape (the meshes of the
        # boolean shapes are the meshes of their operands)
        s isa Union{BooleanUnion,BooleanSubtraction,BooleanIntersection,MultiUnion} && continue
        @test signed_volume(P, F) ≈ capacity(s) rtol = 0.02
    end
    @testset "transformation" begin
        s = Tube{T}(1, 2, 3, 0, 2π)
        t = Transformation3D{T}(10, 0, 0, RotXYZ{T}(π/2, 0, 0))
        m = shape_to_mesh(s; transformation=t)
        P = GeometryBasics.coordinates(m)
        @test signed_volume(P, GeometryBasics.faces(m)) ≈ capacity(s) rtol = 0.02
        @test sum(P) / length(P) ≈ GeometryBasics.Point3{T}(10, 0, 0) atol = 1e-3
        # normals rotated with the shape: along the axis (±y) on the end caps, perpendicular to it
        # on the cylinders
        @test count(n -> abs(n[2]) ≈ 1, m.normal) > 0
        @test all(n -> abs(n[2]) ≈ 1 || abs(n[2]) < 1e-3, m.normal)
    end
end
