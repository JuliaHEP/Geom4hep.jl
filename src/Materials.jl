#---Materials (as G4Isotope, G4Element and G4Material of Geant4)--------------------------------------
#   The quantities are in the units of Geom4hep.Units: lengths in mm, densities in g/cm³ (masses in
#   mg), molar masses in mg/mole (`A = 55.845g/mole`), energies in MeV, temperatures in kelvin.
using Printf: @sprintf
using .Units: g, mole, cm2, cm3, eV, kelvin, atmosphere, Avogadro, amu, fine_structure_const,
              classic_electr_radius, NTP_Temperature, STP_Pressure, universe_mean_density

include("NistData.jl")

#---Isotope-----------------------------------------------------------------------------------------
"""
    Isotope{T}(name, Z, N, A=0; m=0)

Isotope with atomic number `Z`, number of nucleons `N`, molar mass `A` (e.g. `A=55.935g/mole`) and
isomer level `m` (as G4Isotope). When `A` is 0, the mass of the isotope is taken from the NIST data
(natural isotopes only).
"""
struct Isotope{T<:AbstractFloat} <: AbstractMaterial{T}
    name::String
    Z::Int32
    N::Int32
    A::T
    m::Int32
end
function Isotope{T}(name::AbstractString, Z::Integer, N::Integer, A::Real=0; m::Integer=0) where T<:AbstractFloat
    Z >= 1 || throw(ArgumentError("Isotope $name: Z = $Z < 1"))
    N >= Z || throw(ArgumentError("Isotope $name: N = $N < Z = $Z"))
    if A <= 0
        A = _nist_isotope_mass(Z, N)
        A > 0 || throw(ArgumentError("Isotope $name: no NIST mass for Z = $Z, N = $N; give A"))
    end
    Isotope{T}(String(name), Int32(Z), Int32(N), T(A), Int32(m))
end

#---Element-----------------------------------------------------------------------------------------
"""
    Element{T}(name, symbol, Z, A)
    Element{T}(name, symbol, isotopes, abundances)

Chemical element (as G4Element), given by its effective atomic number `Z` and molar mass `A` (e.g.
`Element{T}("Iron", "Fe", 26, 55.845g/mole)`), or by its isotopes and their relative abundances in
number of atoms. In the first case, the isotopes are the natural ones of the NIST data.

The element holds the derived quantities of Geant4: the effective number of nucleons `N`
(`A/(g/mole)`), the Coulomb correction `fcoulomb` and the Tsai factor `radtsai` of the radiation
length (per atom, in mm²), and the mean excitation energy `I` of the NIST data. See also
[`nist_element`](@ref).
"""
struct Element{T<:AbstractFloat} <: AbstractMaterial{T}
    name::String
    symbol::String
    Z::T
    N::T
    A::T
    isotopes::Vector{Isotope{T}}
    abundances::Vector{T}
    I::T
    fcoulomb::T
    radtsai::T
end

function Element{T}(name::AbstractString, symbol::AbstractString, Z::Real, A::Real) where T<:AbstractFloat
    iz = round(Int, Z)
    iz >= 1 || throw(ArgumentError("Element $name: Z = $Z < 1"))
    N = max(A / (g / mole), 1.0)
    N >= Z || throw(ArgumentError("Element $name: N = $N < Z = $Z is not allowed"))
    isotopes, abundances = Isotope{T}[], T[]
    if iz <= length(_NIST_ELEMENTS)
        nist = _NIST_ELEMENTS[iz]
        isempty(symbol) && (symbol = nist.symbol)
        for (n, a, w) in nist.isotopes
            push!(isotopes, Isotope{T}(string(symbol, n), iz, n, a * g / mole))
            push!(abundances, w)
        end
    end
    _element(T, name, symbol, Z, N, A, isotopes, abundances)
end

function Element{T}(name::AbstractString, symbol::AbstractString, isotopes::AbstractVector{<:Isotope},
                    abundances::AbstractVector{<:Real}) where T<:AbstractFloat
    isempty(isotopes) && throw(ArgumentError("Element $name: no isotopes"))
    length(isotopes) == length(abundances) || throw(ArgumentError("Element $name: $(length(isotopes)) isotopes and $(length(abundances)) abundances"))
    Z = first(isotopes).Z
    all(i -> i.Z == Z, isotopes) || throw(ArgumentError("Element $name: isotopes with different Z"))
    w = sum(abundances)
    w > 0 || throw(ArgumentError("Element $name: the sum of the abundances is not positive"))
    ab = T.(abundances ./ w)
    A = sum(ab[i] * isotopes[i].A for i in eachindex(isotopes))
    _element(T, name, symbol, Z, A / (g / mole), A, collect(Isotope{T}, isotopes), ab)
end

function _element(::Type{T}, name, symbol, Z, N, A, isotopes, abundances) where T
    iz = round(Int, Z)
    I = iz <= length(_NIST_ELEMENTS) ? _NIST_ELEMENTS[iz].I * eV : 0.0
    # Coulomb correction factor and Tsai's radiation length factor (Phys Rev. D50 3-1 (1994) 1254)
    az2 = (fine_structure_const * Z)^2
    az4 = az2 * az2
    fcoulomb = (0.0083 * az4 + 0.20206 + 1 / (1 + az2)) * az2 - (0.0020 * az4 + 0.0369) * az4
    if iz <= 4
        Lrad, Lprad = (5.31, 4.79, 4.74, 4.71)[iz], (6.144, 5.621, 5.805, 5.924)[iz]
    else
        logZ3 = log(Z) / 3
        Lrad, Lprad = log(184.15) - logZ3, log(1194.0) - 2logZ3
    end
    radtsai = 4 * fine_structure_const * classic_electr_radius^2 * Z * (Z * (Lrad - fcoulomb) + Lprad)
    Element{T}(String(name), String(symbol), T(Z), T(N), T(A), isotopes, abundances, T(I), T(fcoulomb), T(radtsai))
end

#---Material----------------------------------------------------------------------------------------
"""
    Material{T}(name, Z, A, density; state=:undefined, temperature=NTP_Temperature, pressure=STP_Pressure, I=0)
    Material{T}(name, density; atoms=[...] or fractions=[...], state=:undefined, temperature=NTP_Temperature, pressure=STP_Pressure, I=0, formula="")
    Material{T}(name; density=0, state=:undefined, temperature=NTP_Temperature, pressure=STP_Pressure)

Material of a volume (as G4Material), with the quantities in the units of [`Geom4hep.Units`](@ref)
(the density in g/cm³, e.g. `7.874g/cm3` or `7.874`).

- A material of a single element of atomic number `Z` and molar mass `A` (the element is the NIST
  one for `Z`, as in Geant4).
- A compound given by the number of atoms of each element in a molecule, `atoms=[element => n,
  ...]` (as `G4Material::AddElementByNumberOfAtoms`), e.g. `atoms=[H => 2, O => 1]`.
- A mixture given by mass fractions of elements and materials, `fractions=[element_or_material =>
  w, ...]` (as `AddElementByMassFraction` and `AddMaterial`), e.g. `fractions=[N => 0.755, O =>
  0.232, Ar => 0.013]` or `fractions=[water => 0.5, lead => 0.5]`. The fractions are normalized,
  with a warning when their sum differs from 1 by more than 0.1%.
- A material with only a density (no composition), e.g. for the navigation or `xray` with the mass
  thickness: its radiation and interaction lengths are `NaN` (`Inf` if the density is 0).

`state` is `:solid`, `:liquid`, `:gas` or `:undefined` (solid above 10 mg/cm³, gas below). `I` is
the mean excitation energy; when it is 0, it is computed from the ones of the elements (as
G4IonisParamMat). The derived quantities are computed as in Geant4: see
[`radiation_length`](@ref), [`nuclear_interaction_length`](@ref), [`electron_density`](@ref),
[`mean_excitation_energy`](@ref), [`mass_fractions`](@ref) and [`atoms_per_volume`](@ref). See also
[`nist_material`](@ref).
"""
struct Material{T<:AbstractFloat} <: AbstractMaterial{T}
    name::String
    formula::String
    state::Symbol
    density::T
    temperature::T
    pressure::T
    elements::Vector{Element{T}}
    massfractions::Vector{T}
    natoms::Vector{Int32}             # number of atoms of each element in a molecule (rounded)
    atomspervolume::Vector{T}         # per element
    totatomspervolume::T
    electrondensity::T
    radlen::T
    nuclen::T
    I::T
end

const _kGasThreshold = 10 * 1e-3 * g / cm3          # 10 mg/cm3

function _state(state, density)
    s = Symbol(state)
    s in (:solid, :liquid, :gas, :undefined) || throw(ArgumentError("unknown state $state"))
    s == :undefined ? (density > _kGasThreshold ? :solid : :gas) : s
end

# Material with only a density
function Material{T}(name::AbstractString; density::Real=0, state=:undefined, temperature::Real=NTP_Temperature,
                     pressure::Real=STP_Pressure, I::Real=0, formula::AbstractString="") where T<:AbstractFloat
    unknown = density == 0 ? T(Inf) : T(NaN)
    Material{T}(String(name), String(formula), _state(state, density), T(density), T(temperature), T(pressure),
                Element{T}[], T[], Int32[], T[], zero(T), zero(T), unknown, unknown, T(I))
end

# Material of a single element
function Material{T}(name::AbstractString, Z::Real, A::Real, density::Real; state=:undefined,
                     temperature::Real=NTP_Temperature, pressure::Real=STP_Pressure, I::Real=0,
                     formula::AbstractString="") where T<:AbstractFloat
    density = max(density, universe_mean_density)
    iz = round(Int, Z)
    el = 1 <= iz <= length(_NIST_ELEMENTS) && abs(Z - iz) < 1e-6 ? nist_element(T, iz) :
                                                                   Element{T}("ELM_" * name, name, Z, A)
    _material(T, name, formula, state, density, temperature, pressure, I, [el], [one(T)])
end

# Material made of elements (number of atoms) or of elements and materials (mass fractions)
function Material{T}(name::AbstractString, density::Real; atoms=nothing, fractions=nothing, state=:undefined,
                     temperature::Real=NTP_Temperature, pressure::Real=STP_Pressure, I::Real=0,
                     formula::AbstractString="") where T<:AbstractFloat
    isnothing(atoms) == isnothing(fractions) && throw(ArgumentError("Material $name: give either atoms or fractions"))
    components = something(atoms, fractions)
    isempty(components) && throw(ArgumentError("Material $name: no components"))
    density = max(density, universe_mean_density)
    elements, weights = Element{T}[], T[]
    function add!(el, w)
        k = findfirst(e -> e === el, elements)
        isnothing(k) ? (push!(elements, el); push!(weights, w)) : (weights[k] += w)
    end
    if !isnothing(atoms)                                    # number of atoms in the molecule
        for (el, n) in atoms
            el isa Element || throw(ArgumentError("Material $name: $(el) is not an Element"))
            n isa Real && n > 0 && isinteger(n) || throw(ArgumentError("Material $name: $n atoms of $(el.name)"))
            add!(el, n * el.A)
        end
    else                                                    # mass fractions
        for (c, w) in fractions
            if c isa Element
                add!(c, w)
            elseif c isa Material
                isempty(c.elements) && throw(ArgumentError("Material $name: component $(c.name) has no composition"))
                for (el, f) in zip(c.elements, c.massfractions)
                    add!(el, w * f)
                end
            else
                throw(ArgumentError("Material $name: component $(c) is not an Element or a Material"))
            end
        end
        s = sum(weights)
        abs(1 - s) > 1e-3 && @warn "Material $name: the sum of the mass fractions is $s, not 1"
    end
    _material(T, name, formula, state, density, temperature, pressure, I, elements, weights ./ sum(weights))
end

function _material(::Type{T}, name, formula, state, density, temperature, pressure, I, elements, massfractions) where T
    # molar mass of the "molecule" and number of atoms of each element in it
    Amol = sum(massfractions[i] * elements[i].A for i in eachindex(elements))
    natoms = [Int32(round(massfractions[i] * Amol / elements[i].A)) for i in eachindex(elements)]
    # atoms and electrons per volume
    nv = [Avogadro * density * massfractions[i] / elements[i].A for i in eachindex(elements)]
    totatoms = sum(nv)
    electrons = sum(nv[i] * elements[i].Z for i in eachindex(elements))
    # radiation length
    radinv = sum(nv[i] * elements[i].radtsai for i in eachindex(elements))
    radlen = radinv > 0 ? 1 / radinv : Inf
    # nuclear interaction length (lambda0 = 35 g/cm2 for A=1)
    nilinv = sum(nv[i] * (round(Int, elements[i].Z) == 1 ? elements[i].N : elements[i].N^(2/3)) for i in eachindex(elements))
    nilinv *= amu * cm2 / (35g)
    nuclen = nilinv > 0 ? 1 / nilinv : Inf
    # mean excitation energy: given (NIST DB, GDML <MEE>, user), or from the elements
    if I <= 0 && electrons > 0
        I = exp(sum(nv[i] * elements[i].Z * log(elements[i].I) for i in eachindex(elements)) / electrons)
    end
    Material{T}(String(name), String(formula), _state(state, density), T(density), T(temperature), T(pressure),
                elements, T.(massfractions), natoms, T.(nv), T(totatoms), T(electrons), T(radlen), T(nuclen), T(I))
end

#---Accessors---------------------------------------------------------------------------------------
"Radiation length of a [`Material`](@ref) (mm), with Tsai's formula as in Geant4."
radiation_length(m::Material) = m.radlen
"Nuclear interaction length of a [`Material`](@ref) (mm), with the approximation of Geant4 (35 g/cm² × A^(2/3) per atom)."
nuclear_interaction_length(m::Material) = m.nuclen
"Number of electrons per volume (1/mm³) of a [`Material`](@ref)."
electron_density(m::Material) = m.electrondensity
"Mean excitation energy (MeV) of a [`Material`](@ref) or an [`Element`](@ref)."
mean_excitation_energy(m::Union{Material,Element}) = m.I
"Mass fractions of the elements of a [`Material`](@ref)."
mass_fractions(m::Material) = m.massfractions
"Number of atoms per volume (1/mm³) of each element of a [`Material`](@ref)."
atoms_per_volume(m::Material) = m.atomspervolume

#---NIST elements and materials---------------------------------------------------------------------
const _nist_cache = Dict{Tuple{DataType,Any},Any}()
const _nist_lock = ReentrantLock()

function _nist_isotope_mass(Z, N)
    1 <= Z <= length(_NIST_ELEMENTS) || return 0.0
    for (n, a, _) in _NIST_ELEMENTS[Z].isotopes
        n == N && return a * g / mole
    end
    return 0.0
end

"""
    nist_element([T=Float64,] symbol_or_Z)

Element of the NIST data with its natural isotopic composition (as G4NistManager::FindOrBuildElement),
given by its symbol (e.g. `"Fe"`) or atomic number. The elements are built once and shared.
"""
nist_element(id) = nist_element(Float64, id)
function nist_element(::Type{T}, symbol::AbstractString) where T<:AbstractFloat
    Z = findfirst(e -> e.symbol == symbol, _NIST_ELEMENTS)
    isnothing(Z) && throw(ArgumentError("unknown NIST element $symbol"))
    nist_element(T, Z)
end
function nist_element(::Type{T}, Z::Integer) where T<:AbstractFloat
    1 <= Z <= length(_NIST_ELEMENTS) || throw(ArgumentError("no NIST element with Z = $Z"))
    lock(_nist_lock) do
        get!(_nist_cache, (T, Z)) do
            e = _NIST_ELEMENTS[Z]
            isos = [Isotope{T}(string(e.symbol, n), Z, n, a * g / mole) for (n, a, _) in e.isotopes]
            Element{T}(e.symbol, e.symbol, isos, [w for (_, _, w) in e.isotopes])
        end::Element{T}
    end
end

"""
    nist_material([T=Float64,] name)

Material of the NIST data of Geant4 (as G4NistManager::FindOrBuildMaterial), e.g. `"G4_WATER"`,
`"G4_AIR"`, `"G4_Si"` or `"G4_PbWO4"`: the simple materials `G4_<symbol>` of the elements, and a
selection of compounds and HEP materials (see [`nist_material_names`](@ref)). The materials are
built once and shared.
"""
nist_material(name::AbstractString) = nist_material(Float64, name)
function nist_material(::Type{T}, name::AbstractString) where T<:AbstractFloat
    haskey(_NIST_MATERIALS, name) || throw(ArgumentError("unknown NIST material $name"))
    lock(_nist_lock) do
        get!(_nist_cache, (T, String(name))) do
            d = _NIST_MATERIALS[name]
            comps = [nist_element(T, z) => w for (z, w) in d.components]
            Material{T}(name, d.density; (d.bycount ? (atoms=comps,) : (fractions=comps,))..., state=d.state,
                        temperature=d.temperature, pressure=d.pressure * atmosphere, I=d.I * eV, formula=d.formula)
        end::Material{T}
    end
end

"Names of the materials available with [`nist_material`](@ref)."
nist_material_names() = sort!(collect(keys(_NIST_MATERIALS)))

#---Printing----------------------------------------------------------------------------------------
_fmt(x) = @sprintf("%.4g", x)
function _length(x)
    isfinite(x) || return string(x)
    x >= 1000 ? _fmt(x / 1000) * " m" : x >= 10 ? _fmt(x / 10) * " cm" : _fmt(x) * " mm"
end

Base.show(io::IO, i::Isotope) = print(io, "Isotope(", i.name, ", Z=", i.Z, ", N=", i.N, ", A=", _fmt(i.A / (g / mole)), " g/mole)")
Base.show(io::IO, e::Element) = print(io, "Element(", e.name, " (", e.symbol, "), Z=", _fmt(e.Z), ", A=", _fmt(e.A / (g / mole)), " g/mole)")
Base.show(io::IO, m::Material) = print(io, "Material(", m.name, ", ", _fmt(m.density / (g / cm3)), " g/cm³)")

function Base.show(io::IO, ::MIME"text/plain", e::Element{T}) where T
    println(io, "Element{$T} ", e.name, " (", e.symbol, ")  Z = ", _fmt(e.Z), "  N = ", _fmt(e.N),
            "  A = ", _fmt(e.A / (g / mole)), " g/mole  I = ", _fmt(e.I / eV), " eV")
    for (i, w) in zip(e.isotopes, e.abundances)
        println(io, "  ", rpad(i.name, 6), " Z = ", lpad(i.Z, 3), "  N = ", lpad(i.N, 3), "  A = ",
                lpad(_fmt(i.A / (g / mole)), 7), " g/mole  abundance = ", _fmt(100w), " %")
    end
end

function Base.show(io::IO, ::MIME"text/plain", m::Material{T}) where T
    print(io, "Material{$T} ", m.name)
    isempty(m.formula) || print(io, " (", m.formula, ")")
    println(io, "  ", m.state, "  density = ", _fmt(m.density / (g / cm3)), " g/cm³  T = ",
            _fmt(m.temperature / kelvin), " K  P = ", _fmt(m.pressure / atmosphere), " atm")
    print(io, "  X0 = ", _length(m.radlen), "  λI = ", _length(m.nuclen))
    m.I > 0 && print(io, "  I = ", _fmt(m.I / eV), " eV")
    m.electrondensity > 0 && print(io, "  electrons = ", _fmt(m.electrondensity * cm3), " /cm³")
    for (k, e) in enumerate(m.elements)
        print(io, "\n  ", rpad(e.name, 6), " Z = ", lpad(_fmt(e.Z), 4), "  A = ", lpad(_fmt(e.A / (g / mole)), 7),
              " g/mole  mass fraction = ", lpad(_fmt(100m.massfractions[k]), 6), " %  atoms = ",
              _fmt(m.atomspervolume[k] * cm3), " /cm³")
    end
end
