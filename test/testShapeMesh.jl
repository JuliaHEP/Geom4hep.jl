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
        # a closed surface oriented outwards encloses the volume of the shape (the mesh of a
        # MultiUnion is the mesh of its components)
        s isa MultiUnion && continue
        @test signed_volume(P, F) ≈ capacity(s) rtol = 0.02
    end
    @testset "csg_mesh $(shapename(s))" for s in test_shapes(T)
        s isa Union{BooleanUnion,BooleanSubtraction,BooleanIntersection,MultiUnion} || continue
        m = csg_mesh(s)
        P, F = GeometryBasics.coordinates(m), GeometryBasics.faces(m)
        # closed (no T-junctions): every edge is shared by two triangles, in opposite directions
        edges = Dict{Tuple{Int,Int},Int}()
        for f in F, k in 1:3
            e = (Int(f[k]), Int(f[mod1(k + 1, 3)]))
            edges[e] = get(edges, e, 0) + 1
        end
        @test all(==(1), values(edges))
        @test all(e -> haskey(edges, (e[2], e[1])), keys(edges))
        @test signed_volume(P, F) ≈ capacity(s) rtol = 0.02
    end
    @testset "csg_mesh coplanar faces" begin
        a = Box{T}(1, 1, 1)
        # touching boxes: the common face (x = 1) disappears, and the fragments of the other faces
        # are merged: a box of 6 faces (12 triangles)
        m = csg_mesh(BooleanUnion(a, a, Transformation3D{T}(2, 0, 0)))
        P, F = GeometryBasics.coordinates(m), GeometryBasics.faces(m)
        @test !any(f -> all(i -> abs(P[i][1] - 1) < 1e-5, f), F)
        @test length(F) == 12
        @test signed_volume(P, F) ≈ 16
        # a slot through the box, flush with two of its faces
        m = csg_mesh(BooleanSubtraction(a, Box{T}(0.5, 0.5, 1)))
        @test signed_volume(GeometryBasics.coordinates(m), GeometryBasics.faces(m)) ≈ 6
        # a box with a face in common with the subtracted one: empty
        @test isempty(GeometryBasics.faces(csg_mesh(BooleanSubtraction(a, a))))
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
