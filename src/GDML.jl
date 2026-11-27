using LightXML, Printf
const name = LightXML.name

#--- Unit lookup used when reading GDML (avoids eval of attribute strings)
using .Units: mm, cm, m, km, um, nm, rad, radian, mrad, deg, degree
const _length_units = Dict{String,Float64}("mm" => mm, "cm" => cm, "m" => m, "km" => km,
                                           "um" => um, "micrometer" => um, "nm" => nm)
const _angle_units  = Dict{String,Float64}("rad" => rad, "radian" => radian, "mrad" => mrad,
                                           "deg" => deg, "degree" => degree)

function _parse_unit(table::Dict{String,Float64}, unit::AbstractString)
    haskey(table, unit) || throw(ArgumentError("Unknown or unsupported unit \"$unit\""))
    table[unit]
end

"""
    _unit(attrs, key, default)

Unit factor for the GDML attribute `key` (e.g. "lunit", "aunit", "unit"). Missing attributes take
the GDML defaults (mm for lengths, rad for angles).
"""
function _unit(attrs::AbstractDict, key::AbstractString, default::AbstractString)
    u = get(attrs, key, default)
    haskey(_length_units, u) && return _length_units[u]
    _parse_unit(_angle_units, u)
end
_lunit(attrs) = _unit(attrs, "lunit", "mm")
_aunit(attrs) = _unit(attrs, "aunit", "rad")

# Units of the quantities of the materials (densities, molar masses, energies, temperatures, pressures)
using .Units: g, mg, kg, mole, cm3, m3, eV, keV, MeV, GeV, kelvin, pascal, bar, atmosphere
const _material_units = Dict{String,Float64}(
    "g/cm3" => g/cm3, "mg/cm3" => mg/cm3, "kg/m3" => kg/m3, "g/cm^3" => g/cm3, "kg/m^3" => kg/m3,
    "g/mole" => g/mole, "g/mol" => g/mole, "kg/mole" => kg/mole, "mg/mole" => mg/mole,
    "eV" => eV, "keV" => keV, "MeV" => MeV, "GeV" => GeV,
    "K" => kelvin, "kelvin" => kelvin,
    "pascal" => pascal, "hep_pascal" => pascal, "Pa" => pascal, "bar" => bar, "atmosphere" => atmosphere, "atm" => atmosphere)

# Value of a quantity of a material: an element <X value= unit=/> (with the default unit of the
# GDML schema), or a reference <Xref ref=/> to a <quantity> of the <define> section
function _material_quantity(dicts, e::XMLElement, default::AbstractString)
    attrs = attributes_dict(e)
    if haskey(attrs, "ref")
        haskey(dicts.quantities, attrs["ref"]) || throw(ArgumentError("GDML: quantity $(attrs["ref"]) is not defined"))
        return dicts.quantities[attrs["ref"]]
    end
    parse(Float64, attrs["value"]) * _parse_unit(_material_units, get(attrs, "unit", default))
end

struct GDMLDicts{T<:AbstractFloat}
    materials::Dict{String,AbstractMaterial{T}}
    solids::Dict{String,AbstractShape{T}}
    volumes::Dict{String,Volume{T}}
    positions::Dict{String, Vector3{T}}
    rotations::Dict{String, RotXYZ{T}}
    quantities::Dict{String, Float64}
    function GDMLDicts{T}() where T<:AbstractFloat
        new(Dict{String,AbstractMaterial{T}}(),
            Dict{String,AbstractShape{T}}(),
            Dict{String,Volume{T}}(),
            Dict{String, Vector3{T}}(),
            Dict{String, RotXYZ{T}}(),
            Dict{String, Float64}()
        )
    end
end

# process <define/>
function fill_define!(dicts::GDMLDicts{T}, element::XMLElement) where T<:AbstractFloat
    (; positions, rotations) = dicts
    for c in child_nodes(element)
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            attrs = attributes_dict(e)
            if elemname == "position"
                unit = _unit(attrs, "unit", "mm")
                positions[attrs["name"]] = Vector3{T}(parse(T, attrs["x"]) * unit, 
                                                      parse(T, attrs["y"]) * unit, 
                                                      parse(T, attrs["z"]) * unit)
            elseif elemname == "rotation"
                unit = _unit(attrs, "unit", "rad")
                rotations[attrs["name"]] = RotZYX{T}(parse(T, attrs["z"]) * unit, 
                                                     parse(T, attrs["y"]) * unit, 
                                                     parse(T, attrs["x"]) * unit)
            elseif elemname == "quantity"
                u = get(attrs, "unit", "")
                f = haskey(_material_units, u) ? _material_units[u] : haskey(_length_units, u) ? _length_units[u] :
                    haskey(_angle_units, u) ? _angle_units[u] : 1.0
                dicts.quantities[attrs["name"]] = parse(Float64, attrs["value"]) * f
            elseif elemname in ("matrix", "constant", "variable")
                continue                         # e.g. tables of optical properties: not used
            else
                @warn "GDML: define $elemname not yet supported" maxlog=1
            end
        end
    end
end

# process <materials/> (as G4GDMLReadMaterials)
function fill_materials!(dicts::GDMLDicts{T}, element::XMLElement) where T<:AbstractFloat
    (; materials) = dicts
    children(e) = (XMLElement(c) for c in child_nodes(e) if is_elementnode(c))
    for e in children(element)
        elemname = name(e)
        attrs = attributes_dict(e)
        matname = attrs["name"]
        if elemname == "isotope"
            A = 0.0
            for c in children(e)
                name(c) in ("atom", "atomref") && (A = _material_quantity(dicts, c, "g/mole"))
            end
            materials[matname] = Isotope{T}(matname, parse(Int, attrs["Z"]), parse(Int, attrs["N"]), A)
        elseif elemname == "element"
            symbol = get(attrs, "formula", "")
            A = 0.0
            isotopes, abundances = Isotope{T}[], Float64[]
            for c in children(e)
                if name(c) in ("atom", "atomref")
                    A = _material_quantity(dicts, c, "g/mole")
                elseif name(c) == "fraction"
                    aa = attributes_dict(c)
                    push!(isotopes, materials[aa["ref"]])
                    push!(abundances, parse(Float64, aa["n"]))
                end
            end
            if !isempty(isotopes)
                materials[matname] = Element{T}(matname, symbol, isotopes, abundances)
            else
                Z = parse(Float64, attrs["Z"])
                A > 0 || (A = nist_element(T, round(Int, Z)).A)       # natural element of the NIST data
                materials[matname] = Element{T}(matname, symbol, Z, A)
            end
        elseif elemname == "material"
            kw = Dict{Symbol,Any}(:state => Symbol(get(attrs, "state", "undefined")),
                                  :formula => get(attrs, "formula", ""))
            density, A = 0.0, 0.0
            fractions, atoms = Pair[], Pair[]
            for c in children(e)
                cname = name(c)
                if cname in ("D", "Dref")
                    density = _material_quantity(dicts, c, "g/cm3")
                elseif cname in ("T", "Tref")
                    kw[:temperature] = _material_quantity(dicts, c, "K")
                elseif cname in ("P", "Pref")
                    kw[:pressure] = _material_quantity(dicts, c, "pascal")
                elseif cname in ("MEE", "MEEref")
                    kw[:I] = _material_quantity(dicts, c, "eV")
                elseif cname in ("atom", "atomref")
                    A = _material_quantity(dicts, c, "g/mole")
                elseif cname == "fraction"
                    aa = attributes_dict(c)
                    push!(fractions, _gdml_component(dicts, aa["ref"]) => parse(Float64, aa["n"]))
                elseif cname == "composite"
                    aa = attributes_dict(c)
                    push!(atoms, _gdml_component(dicts, aa["ref"]) => parse(Int, aa["n"]))
                end
            end
            Z = parse(Float64, get(attrs, "Z", "0"))
            materials[matname] = if !isempty(atoms) && !isempty(fractions)
                throw(ArgumentError("GDML: material $matname has both <composite> and <fraction>"))
            elseif !isempty(atoms)
                Material{T}(matname, density; atoms, kw...)
            elseif !isempty(fractions)
                Material{T}(matname, density; fractions, kw...)
            elseif Z >= 1
                Material{T}(matname, Z, A, density; kw...)
            else                                       # e.g. a dummy material with Z=0 (DD4hep)
                Material{T}(matname; density, kw...)
            end
        end
    end
end

# Element or material referenced by a GDML material: defined in the file, or from the NIST data
# (e.g. "G4_AIR" or "Fe", as G4GDMLReadMaterials)
function _gdml_component(dicts::GDMLDicts{T}, ref::AbstractString) where T
    haskey(dicts.materials, ref) && return dicts.materials[ref]
    haskey(_NIST_MATERIALS, ref) && return (dicts.materials[ref] = nist_material(T, ref))
    any(e -> e.symbol == ref, _NIST_ELEMENTS) && return (dicts.materials[ref] = nist_element(T, ref))
    throw(ArgumentError("GDML: material or element $ref is not defined"))
end

function _gdml_material(dicts::GDMLDicts{T}, ref::AbstractString) where T
    m = _gdml_component(dicts, ref)
    m isa Material || throw(ArgumentError("GDML: $ref is not a material"))
    return m
end

# process <solids/>
function fill_solids!(dicts::GDMLDicts{T}, element::XMLElement) where T<:AbstractFloat
    (; solids) = dicts
    unsupported = Dict{String,Int}()
    for c in child_nodes(element)
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            attrs = attributes_dict(e)
            if elemname == "box"
                lunit = _lunit(attrs)
                solids[attrs["name"]] = Box{T}(parse(T, attrs["x"]) * lunit / 2, 
                                               parse(T, attrs["y"]) * lunit / 2, 
                                               parse(T, attrs["z"]) * lunit / 2)
            elseif elemname == "trd"
                lunit = _lunit(attrs)
                solids[attrs["name"]] = Trd{T}(parse(T, attrs["x1"]) * lunit / 2, 
                                               parse(T, attrs["x2"]) * lunit / 2, 
                                               parse(T, attrs["y1"]) * lunit / 2,
                                               parse(T, attrs["y2"]) * lunit / 2,
                                               parse(T, attrs["z"])  * lunit / 2)
            elseif elemname == "trap"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                solids[attrs["name"]] = Trap{T}(parse(T, attrs["z"]) * lunit / 2, 
                                                parse(T, attrs["theta"]) * aunit, 
                                                parse(T, attrs["phi"]) * aunit,
                                                parse(T, attrs["y1"]) * lunit / 2,
                                                parse(T, attrs["x1"])  * lunit / 2,
                                                parse(T, attrs["x2"])  * lunit / 2,
                                                parse(T, attrs["alpha1"])  * aunit,
                                                parse(T, attrs["y2"])  * lunit / 2,
                                                parse(T, attrs["x3"])  * lunit / 2,
                                                parse(T, attrs["x4"])  * lunit / 2,
                                                parse(T, attrs["alpha2"])  * aunit)
            elseif elemname == "tube"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                solids[attrs["name"]] = Tube{T}(parse(T, attrs["rmin"]) * lunit, 
                                                parse(T, attrs["rmax"]) * lunit, 
                                                parse(T, attrs["z"]) * lunit / 2,
                                                parse(T, attrs["startphi"]) * aunit,
                                                parse(T, attrs["deltaphi"]) * aunit)
            elseif elemname == "sphere"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                attr(k, default) = haskey(attrs, k) ? parse(T, attrs[k]) : T(default)
                solids[attrs["name"]] = Sphere{T}(attr("rmin", 0) * lunit, parse(T, attrs["rmax"]) * lunit,
                                                  attr("startphi", 0) * aunit, attr("deltaphi", 2π/aunit) * aunit,
                                                  attr("starttheta", 0) * aunit, attr("deltatheta", π/aunit) * aunit)
            elseif elemname == "orb"
                solids[attrs["name"]] = Sphere{T}(0, parse(T, attrs["r"]) * _lunit(attrs))
            elseif elemname == "opticalsurface"
                continue                         # optical surfaces are not part of the geometry
            elseif elemname == "cone"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                solids[attrs["name"]] = Cone{T}(parse(T, attrs["rmin1"]) * lunit, 
                                                parse(T, attrs["rmax1"]) * lunit,
                                                parse(T, attrs["rmin2"]) * lunit,
                                                parse(T, attrs["rmax2"]) * lunit, 
                                                parse(T, attrs["z"]) * lunit / 2,
                                                parse(T, attrs["startphi"]) * aunit,
                                                parse(T, attrs["deltaphi"]) * aunit)
            elseif elemname == "polycone"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                rmax, rmin, z = Vector{T}(), Vector{T}(), Vector{T}()
                for cc in child_nodes(e)
                    if name(cc) == "zplane"
                        aa = attributes_dict(XMLElement(cc))
                        push!(rmax, parse(T, aa["rmax"]) * lunit)
                        push!(rmin, parse(T, aa["rmin"]) * lunit)
                        push!(z, parse(T, aa["z"]) * lunit)
                    end
                end
                N = length(rmax)
                solids[attrs["name"]] = Polycone{T}(rmin, rmax, z,
                                                    parse(T, attrs["startphi"]) * aunit,
                                                    parse(T, attrs["deltaphi"]) * aunit)
            elseif elemname == "torus"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                solids[attrs["name"]] = Torus{T}(parse(T, get(attrs, "rmin", "0")) * lunit,
                                                 parse(T, attrs["rmax"]) * lunit,
                                                 parse(T, attrs["rtor"]) * lunit,
                                                 parse(T, get(attrs, "startphi", "0")) * aunit,
                                                 parse(T, get(attrs, "deltaphi", "360")) * aunit)
            elseif elemname == "polyhedra"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                rmax, rmin, z = Vector{T}(), Vector{T}(), Vector{T}()
                for cc in child_nodes(e)
                    if name(cc) == "zplane"
                        aa = attributes_dict(XMLElement(cc))
                        push!(rmax, parse(T, aa["rmax"]) * lunit)
                        push!(rmin, parse(T, get(aa, "rmin", "0")) * lunit)
                        push!(z, parse(T, aa["z"]) * lunit)
                    end
                end
                solids[attrs["name"]] = Polyhedra{T}(parse(Int, attrs["numsides"]), rmin, rmax, z,
                                                     parse(T, attrs["startphi"]) * aunit,
                                                     parse(T, attrs["deltaphi"]) * aunit)
            elseif elemname == "cutTube"
                lunit = _lunit(attrs)
                aunit = _aunit(attrs)
                solids[attrs["name"]] = CutTube{T}(parse(T, attrs["rmin"]) * lunit,
                                                   parse(T, attrs["rmax"]) * lunit,
                                                   parse(T, attrs["z"]) * lunit / 2,   # GDML gives the full length
                                                   parse(T, attrs["startphi"]) * aunit,
                                                   parse(T, attrs["deltaphi"]) * aunit,
                                                   Vector3{T}(parse(T, attrs["lowX"]), parse(T, attrs["lowY"]), parse(T, attrs["lowZ"])),
                                                   Vector3{T}(parse(T, attrs["highX"]), parse(T, attrs["highY"]), parse(T, attrs["highZ"])))
            elseif elemname in ("union", "subtraction", "intersection") 
                first = nothing
                second = nothing
                transformation = read_transformation(dicts, e)
                for cc in child_nodes(e)
                    if is_elementnode(cc)
                        aa = attributes_dict(XMLElement(cc))
                        if name(cc) == "first"
                            first = _get_solid(solids, aa["ref"], attrs["name"])
                        elseif name(cc) == "second"
                            second = _get_solid(solids, aa["ref"], attrs["name"])
                        end
                    end
                end
                if isnothing(first) || isnothing(second)
                    # Broken GDML: keep whatever operand is valid (or nothing at all for intersections)
                    solids[attrs["name"]] = elemname == "intersection" || isnothing(first) ? NoShape{T}() : first
                elseif elemname == "union"
                    solids[attrs["name"]] = BooleanUnion(first, second, transformation)
                elseif elemname == "subtraction"
                    solids[attrs["name"]] = BooleanSubtraction(first, second, transformation)
                elseif elemname == "intersection"
                    solids[attrs["name"]] = BooleanIntersection(first, second, transformation)
                end
            else
                unsupported[elemname] = get(unsupported, elemname, 0) + 1
                solids[attrs["name"]] = NoShape{T}()
            end
        end
    end
    isempty(unsupported) || @warn "GDML: unsupported shapes replaced by NoShape" unsupported
end

function _get_solid(solids::Dict, ref::AbstractString, owner::AbstractString)
    haskey(solids, ref) && return solids[ref]
    @warn "GDML: solid \"$ref\" referenced by boolean \"$owner\" is not defined; ignoring this operand"
    return nothing
end

# process <structure/>
function fill_volumes!(dicts::GDMLDicts{T}, element::XMLElement) where T<:AbstractFloat
    (; materials, solids, volumes, positions, rotations) = dicts
    for c in child_nodes(element) #--- loop over volumes
        if is_elementnode(c)
            e = XMLElement(c)
            elemname = name(e)
            attrs = attributes_dict(e)
            if elemname in ("volume", "assembly")
                volname = attrs["name"]
                shape = nothing
                material = nothing
                pvols = PlacedVolume{T}[]
                for cc in child_nodes(e)  #--- loop over ref to solid, material and daughters
                    if is_elementnode(cc)
                        aa = attributes_dict(XMLElement(cc))
                        if  name(cc) == "solidref"
                            shape = solids[aa["ref"]]
                        elseif name(cc) == "materialref"
                            material = _gdml_material(dicts, aa["ref"])
                        elseif name(cc) == "physvol"
                            daughter = nothing
                            transformation = read_transformation(dicts, XMLElement(cc))
                            for ccc in child_nodes(XMLElement(cc))
                                if is_elementnode(ccc)
                                    aa = attributes_dict(XMLElement(ccc))
                                    if name(ccc) == "volumeref"
                                        daughter = volumes[aa["ref"]]
                                    end
                                end
                            end
                            push!(pvols, PlacedVolume{T}(-1, transformation, daughter))
                        end
                    end
                end
                volume = elemname == "volume" ? Volume{T}(volname, shape, material) : Assembly{T}(volname)
                for pvol in pvols
                    place_daughter!(volume, pvol.transformation, pvol.volume)
                end
                volumes[volname] = volume
            end
        end
    end
end

# fill Transformation3D
function read_transformation(dicts::GDMLDicts{T}, element::XMLElement) where T
    (; positions, rotations) = dicts
    position = Vector3{T}(0,0,0)
    rotation = RotZYX{T}(0,0,0)
    for c in child_nodes(element)
        if is_elementnode(c)
            aa = attributes_dict(XMLElement(c))
            if name(c) == "position"
                unit = _unit(aa, "unit", "mm")
                position = Vector3{T}(parse(T, aa["x"]) * unit, 
                                      parse(T, aa["y"]) * unit,
                                      parse(T, aa["z"]) * unit)
            elseif name(c) == "rotation"
                unit = _unit(aa, "unit", "rad")
                rotation = RotZYX{T}(parse(T, aa["z"]) * unit, 
                                     parse(T, aa["y"]) * unit, 
                                     parse(T, aa["x"]) * unit)
            elseif name(c) == "positionref"
                position = positions[aa["ref"]]
            elseif name(c) == "rotationref"
                rotation = rotations[aa["ref"]]
            end
        end
    end
    Transformation3D{T}(RotMatrix3{T}(rotation), position)
end


# process the full file
"""
    read_gdml(filename, T=Float64) -> Volume{T}

Read a GDML file and return its world volume, with the shapes, materials and volume hierarchy it
defines, with floating point type `T`. Lengths are converted to mm and angles to radians. Solids
that are not supported are replaced by [`NoShape`](@ref) with a warning.
"""
function read_gdml(fname::String, ::Type{T}=Float64) where T
    dicts = GDMLDicts{T}()
    world = nothing
    xdoc = parse_file(fname)
    xroot = root(xdoc) 
    if name(xroot) != "gdml"
        throw(ArgumentError("File $fname is not a GDML file"))
    end
    for c in child_nodes(xroot)  # c is an instance of XMLNode
        if is_elementnode(c)
            e = XMLElement(c)    # this makes an XMLElement instance
            if name(e) == "materials"
                fill_materials!(dicts, e)
            elseif name(e) == "solids"
                fill_solids!(dicts, e)
            elseif name(e) == "structure"
                fill_volumes!(dicts, e)
            elseif name(e) == "define"
                fill_define!(dicts, e)
            elseif name(e) == "setup"
                for cc in child_nodes(e)
                    if is_elementnode(cc)
                        ee = XMLElement(cc)
                        attrs = attributes_dict(ee)
                        if name(ee) == "world"
                           world = dicts.volumes[attrs["ref"]] 
                        end
                    end
                end
            end
        end
    end
    free(xdoc)
    return world
end

