using Geom4hep.Units

@testset "Units" begin
    @test g / cm3 == 1                       # densities in g/cm3
    @test kg / m3 ≈ 1e-3
    @test g / mole == 1000
    @test 1 / Geom4hep.Units.fine_structure_const ≈ 137.036 rtol = 1e-5
    @test Geom4hep.Units.classic_electr_radius ≈ 2.8179e-12 rtol = 1e-4      # mm
    @test Geom4hep.Units.amu / g ≈ 1.66054e-24 rtol = 1e-5
    @test atmosphere / bar ≈ 1.01325
end

@testset "Materials{$T}" for T in (Float64, Float32)
    @testset "Isotopes and elements" begin
        fe56 = Isotope{T}("Fe56", 26, 56)                                   # NIST mass
        @test fe56.A ≈ T(55.9349g / mole) rtol = 1e-5
        @test_throws ArgumentError Isotope{T}("X", 26, 30)                  # not a natural isotope
        @test Isotope{T}("X", 26, 30, 30g / mole).A ≈ T(30g / mole)
        @test_throws ArgumentError Isotope{T}("X", 26, 20)                  # N < Z

        fe = nist_element(T, "Fe")
        @test fe === nist_element(T, 26)                                    # built once
        @test fe.Z == 26
        @test fe.A ≈ T(55.845g / mole) rtol = 1e-4
        @test fe.N ≈ fe.A / T(g / mole)
        @test length(fe.isotopes) == 4 && sum(fe.abundances) ≈ 1
        @test mean_excitation_energy(fe) ≈ T(286eV)

        # by Z and A, with the natural isotopes
        o = Element{T}("Oxygen", "O", 8, 16.00g / mole)
        @test o.Z == 8 && o.A ≈ T(16g / mole) && length(o.isotopes) == 3
        # from isotopes (abundances normalized)
        n = Element{T}("N", "N", [Isotope{T}("N14", 7, 14, 14.0031g / mole), Isotope{T}("N15", 7, 15, 15.0001g / mole)], [99.632, 0.368])
        @test n.A ≈ T((0.99632 * 14.0031 + 0.00368 * 15.0001)g / mole)
        @test sum(n.abundances) ≈ 1
        @test_throws ArgumentError Element{T}("X", "X", [fe56, Isotope{T}("O16", 8, 16)], [1, 1])
        @test_throws ArgumentError Element{T}("X", "X", 0, 1g / mole)
    end

    @testset "Constructors" begin
        H, O, C = nist_element(T, "H"), nist_element(T, "O"), nist_element(T, "C")
        # by number of atoms
        water = Material{T}("Water", 1.0g / cm3; atoms=[H => 2, O => 1], I=78eV)
        @test water.natoms == [2, 1]
        @test mass_fractions(water) ≈ T[2H.A, O.A] / (2H.A + O.A)
        @test water.state == :solid                                         # undefined, dense
        @test electron_density(water) ≈ T(3.343e23 / cm3) rtol = 1e-3
        # by mass fractions, with the same element twice
        m = Material{T}("M", 1.0; fractions=[H => 0.1, O => 0.6, H => 0.3], state=:liquid)
        @test length(m.elements) == 2 && mass_fractions(m) ≈ T[0.4, 0.6]
        @test m.state == :liquid
        # mixture of materials
        mix = Material{T}("Mix", 1.0; fractions=[water => 0.5, Material{T}("Carbon", 6, 12.011g / mole, 2.0) => 0.5])
        @test length(mix.elements) == 3
        @test sum(mass_fractions(mix)) ≈ 1
        @test mass_fractions(mix)[3] ≈ T(0.5)                               # carbon
        # single element: the NIST element
        al = Material{T}("Aluminium", 13, 26.98g / mole, 2.699g / cm3)
        @test al.elements[1] === nist_element(T, 13)
        @test radiation_length(al) ≈ T(8.897cm) rtol = 2e-3
        # the fractions are normalized, with a warning
        @test_logs (:warn, r"sum of the mass fractions") Material{T}("W", 1.0; fractions=[H => 0.5, O => 0.6])
        @test_throws ArgumentError Material{T}("X", 1.0; atoms=[H => 2, O => 0.5])   # not a number of atoms
        @test_throws ArgumentError Material{T}("X", 1.0; atoms=[water => 1])           # not an element
        @test_throws ArgumentError Material{T}("X", 1.0; atoms=[H => 2], fractions=[O => 1.0])
        @test_throws ArgumentError Material{T}("X", 1.0)
        @test_throws ArgumentError Material{T}("X", 1.0; fractions=Pair[])
        # an integer mass fraction is still a mass fraction
        ho = @test_logs (:warn, r"sum of the mass fractions") Material{T}("HO", 1.0; fractions=[H => 1, O => 1])
        @test mass_fractions(ho) ≈ T[0.5, 0.5]
        # only a density: unknown lengths (infinite for a vacuum)
        @test isnan(radiation_length(Material{T}("iron"; density=7.87)))
        @test isinf(radiation_length(Material{T}("vacuum"; density=0)))
        @test Material{T}("air"; density=0.0012).state == :gas
    end

    @testset "NIST materials" begin
        # radiation lengths of the PDG, nuclear interaction lengths and densities of electrons of Geant4
        for (name, X0, λ) in (("G4_WATER", 36.08cm, 75.375cm), ("G4_Pb", 0.5612cm, 18.248cm),
                              ("G4_Fe", 1.757cm, NaN), ("G4_Si", 9.370cm, NaN), ("G4_Cu", 1.436cm, NaN),
                              ("G4_AIR", 303.9m, 710.1m), ("G4_lAr", 14.0cm, NaN), ("G4_POLYSTYRENE", 41.31cm, NaN))
            mat = nist_material(T, name)
            @test radiation_length(mat) ≈ T(X0) rtol = 2e-3
            isnan(λ) || @test nuclear_interaction_length(mat) ≈ T(λ) rtol = 1e-3
        end
        @test nist_material(T, "G4_WATER") === nist_material(T, "G4_WATER")
        @test mean_excitation_energy(nist_material(T, "G4_WATER")) ≈ T(78eV)
        @test nist_material(T, "G4_WATER").formula == "H_2O"
        @test nist_material(T, "G4_AIR").state == :gas
        # mean excitation energy from the elements when not in the NIST data (Bragg's rule)
        pbwo4 = nist_material(T, "G4_PbWO4")
        @test 286eV < mean_excitation_energy(pbwo4) / T(1) < 823eV
        galactic = nist_material(T, "G4_Galactic")
        @test galactic.density ≈ T(1e-25) && galactic.temperature ≈ T(2.73)
        @test radiation_length(galactic) > T(1e20)
        @test "G4_PbWO4" in nist_material_names() && "G4_Cf" in nist_material_names()
        @test_throws ArgumentError nist_material(T, "G4_UNKNOWN")
        @test_throws ArgumentError nist_element(T, "Xx")
    end

    @testset "GDML materials" begin
        gdml = """
        <?xml version="1.0" encoding="UTF-8"?>
        <gdml>
          <define>
            <quantity name="rho" type="density" value="2.0" unit="g/cm3"/>
          </define>
          <materials>
            <element name="Hyd" formula="H" Z="1"> <atom unit="g/mole" value="1.008"/> </element>
            <isotope name="O16" N="16" Z="8"> <atom unit="g/mole" value="15.9949"/> </isotope>
            <element name="Oxy" formula="O"> <fraction n="1.0" ref="O16"/> </element>
            <material name="Water" state="liquid">
              <MEE unit="eV" value="78"/> <D unit="kg/m3" value="1000"/>
              <composite n="2" ref="Hyd"/> <composite n="1" ref="Oxy"/>
            </material>
            <material name="WaterAndLead">
              <Dref ref="rho"/>
              <fraction n="0.5" ref="Water"/> <fraction n="0.5" ref="G4_Pb"/>
            </material>
            <material name="Gold" Z="79"> <D value="19.3"/> <atom value="196.97"/> <T unit="K" value="77"/> </material>
            <material name="Dummy" Z="0"> <D value="0"/> <atom value="0"/> </material>
          </materials>
          <solids>
            <box name="WorldBox" x="100" y="100" z="100"/>
            <box name="B" x="10" y="10" z="10"/>
          </solids>
          <structure>
            <volume name="Inner"> <materialref ref="WaterAndLead"/> <solidref ref="B"/> </volume>
            <volume name="World">
              <materialref ref="G4_AIR"/> <solidref ref="WorldBox"/>
              <physvol> <volumeref ref="Inner"/> </physvol>
            </volume>
          </structure>
          <setup name="Default" version="1.0"> <world ref="World"/> </setup>
        </gdml>
        """
        file = tempname() * ".gdml"
        write(file, gdml)
        world = read_gdml(file, T)
        rm(file)
        @test world.material.name == "G4_AIR"                              # from the NIST data
        mix = world.daughters[1].volume.material
        @test mix.density ≈ T(2.0)
        @test length(mix.elements) == 3                                     # H, O and Pb
        @test sum(mass_fractions(mix)) ≈ 1
        @test any(e -> e.symbol == "Pb", mix.elements)
    end

    @testset "GDML materials (details)" begin
        dicts = Geom4hep.GDMLDicts{T}()
        xml = Geom4hep.LightXML.parse_string("""
          <materials>
            <element name="Hyd" formula="H" Z="1"> <atom unit="g/mole" value="1.008"/> </element>
            <element name="Oxy" formula="O" Z="8"> <atom value="16.0"/> </element>
            <material name="Water" state="liquid">
              <MEE unit="eV" value="78"/> <D unit="kg/m3" value="1000"/>
              <composite n="2" ref="Hyd"/> <composite n="1" ref="Oxy"/>
            </material>
            <material name="Gold" Z="79"> <D value="19.3"/> <atom value="196.97"/> <T unit="K" value="77"/> <P unit="bar" value="2"/> </material>
            <material name="Dummy" Z="0"> <D value="0"/> <atom value="0"/> </material>
          </materials>""")
        Geom4hep.fill_materials!(dicts, Geom4hep.LightXML.root(xml))
        water, gold, dummy = dicts.materials["Water"], dicts.materials["Gold"], dicts.materials["Dummy"]
        @test water.density ≈ 1 && water.state == :liquid && water.natoms == [2, 1]
        @test mean_excitation_energy(water) ≈ T(78eV)
        @test radiation_length(water) ≈ T(36.08cm) rtol = 2e-3
        @test gold.temperature ≈ 77 && gold.pressure ≈ T(2bar) && gold.elements[1].Z == 79
        @test isinf(radiation_length(dummy)) && isempty(dummy.elements)
    end

    @testset "X-ray of the material budget" begin
        world = Volume{T}("world", Box{T}(200, 200, 200), nist_material(T, "G4_Galactic"))
        place_daughter!(world, Transformation3D{T}(0, 0, 0), Volume{T}("slab", Box{T}(100, 100, 50), nist_material(T, "G4_Pb")))
        model = close_geometry(world)
        pb = nist_material(T, "G4_Pb")
        for (q, expected) in ((:mass, 100 * pb.density), (:radlen, 100 / radiation_length(pb)), (:nuclen, 100 / nuclear_interaction_length(pb)))
            image = xray(model, 3, 20; quantity=q)[3]
            @test image[10, 10] ≈ expected rtol = 1e-4
            @test image[1, 1] < 1e-6 * expected                             # only G4_Galactic
        end
        @test_throws ArgumentError xray(model, 3, 4; quantity=:foo)
    end
end
