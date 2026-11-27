#!/usr/bin/env python3
"""
Generate src/NistData.jl from the Geant4 NIST builders.

Usage:
    python3 -I utils/generate_nist_data.py <G4NistElementBuilder.cc> <G4NistMaterialBuilder.cc> > src/NistData.jl

The sources are in geant4/source/materials/src of a Geant4 release (see GEANT4_VERSION). The
elements are taken with their natural isotopic composition (isotopes with a non-zero abundance),
and the materials of MATERIALS (the simple G4_<X> materials of all the elements and a selection of
compounds), with the conventions of G4NistMaterialBuilder (density in g/cm3, mean excitation energy
in eV, 0 when it is not defined).
"""
import math
import re
import sys

GEANT4_VERSION = "11.3.0"

# CLHEP constants of this Geant4 version (CLHEP/Units/PhysicalConstants.h), in MeV
electron_mass_c2 = 0.510998910
proton_mass_c2 = 938.272013
amu_c2 = 931.494028
GeV = 1000.0
eV = 1e-6

# Compounds and HEP materials (the G4_<X> materials of all the elements are always included)
MATERIALS = [
    # NIST compounds
    "G4_AIR", "G4_WATER", "G4_WATER_VAPOR", "G4_ALUMINUM_OXIDE", "G4_BARIUM_FLUORIDE", "G4_BGO",
    "G4_BUTANE", "G4_CALCIUM_FLUORIDE", "G4_CARBON_DIOXIDE", "G4_CESIUM_IODIDE", "G4_CONCRETE",
    "G4_ETHANE", "G4_GLASS_PLATE", "G4_KAPTON", "G4_LEAD_OXIDE", "G4_LITHIUM_FLUORIDE",
    "G4_METHANE", "G4_MYLAR", "G4_NYLON-6-6", "G4_PARAFFIN", "G4_PLASTIC_SC_VINYLTOLUENE",
    "G4_PLEXIGLASS", "G4_POLYCARBONATE", "G4_POLYETHYLENE", "G4_POLYPROPYLENE", "G4_POLYSTYRENE",
    "G4_POLYVINYL_CHLORIDE", "G4_PROPANE", "G4_Pyrex_Glass", "G4_SILICON_DIOXIDE",
    "G4_SODIUM_IODIDE", "G4_TEFLON", "G4_BONE_COMPACT_ICRU", "G4_MUSCLE_STRIATED_ICRU",
    "G4_TISSUE_SOFT_ICRP",
    # HEP and nuclear materials
    "G4_lH2", "G4_lN2", "G4_lO2", "G4_lAr", "G4_lKr", "G4_lXe", "G4_PbWO4", "G4_Galactic",
    "G4_GRAPHITE_POROUS", "G4_BRASS", "G4_BRONZE", "G4_STAINLESS-STEEL",
]

STATES = {"kStateSolid": "solid", "kStateLiquid": "liquid", "kStateGas": "gas",
          "kStateUndefined": "undefined"}


def number(s):
    return float(s.strip().rstrip("fF"))


def parse_elements(src):
    """Elements with their natural isotopes: [(symbol, Z, [(N, A, abundance), ...]), ...]"""
    arrays = {}
    for m in re.finditer(r"G4(?:int|double)\s+(\w+)\[\d+\]\s*=\s*\{([^}]*)\}", src):
        arrays[m.group(1)] = [number(x) for x in m.group(2).split(",") if x.strip()]
    # H and He masses are recomputed in Geant4 to be consistent with its particle masses
    def binding(Z):
        return (14.4381 * Z**2.39 + 1.55468e-6 * Z**5.35) * eV
    arrays["HA"][0] = (proton_mass_c2 + electron_mass_c2 - binding(1)) / amu_c2
    arrays["HA"][1] = (1.875613 * GeV + electron_mass_c2 - binding(1)) / amu_c2
    arrays["HA"][2] = (2.80925 * GeV + electron_mass_c2 - binding(1)) / amu_c2
    arrays["HeA"][0] = (2.80923 * GeV + 2.0 * electron_mass_c2 - binding(2)) / amu_c2
    arrays["HeA"][1] = (3.727417 * GeV + 2.0 * electron_mass_c2 - binding(2)) / amu_c2

    elements = []
    for m in re.finditer(r'AddElement\("(\w+)",\s*(\d+),\s*(\d+),\s*\*(\w+),\s*\*(\w+),\s*\*(\w+),\s*\*(\w+)\)', src):
        sym, Z, nc, aN, aA, _, aW = m.groups()
        Z, nc = int(Z), int(nc)
        N0 = int(arrays[aN][0])
        A, W = arrays[aA][:nc], arrays[aW][:nc]
        wsum = sum(W)
        isotopes = [(N0 + i, A[i], W[i] / wsum) for i in range(nc) if W[i] > 0.0]
        elements.append((sym, Z, isotopes))
    return elements


def parse_materials(src, symbols):
    """Materials in the order of G4NistMaterialBuilder: name -> dict"""
    zof = {s: z for s, z in symbols}
    materials = {}
    current = None
    call = re.compile(r'(AddMaterial|AddElementByAtomCount|AddElementByWeightFraction)\(([^;]*)\);'
                      r'|chFormulas\[nMaterials-1\]\s*=\s*"([^"]*)"'
                      r'|AddGas\("(\w+)",\s*([^,]+),\s*([^)]+)\)')
    body = src[src.index("void G4NistMaterialBuilder::NistSimpleMaterials()"):]
    body = re.sub(r"/\*.*?\*/", "", body, flags=re.S)       # C and C++ comments
    body = re.sub(r"//[^\n]*", "", body)
    for m in call.finditer(body):
        if m.group(3) is not None:
            current["formula"] = m.group(3)
            continue
        if m.group(4) is not None:          # AddGas: temperature and pressure of G4_Galactic
            materials[m.group(4)]["gas"] = (m.group(5).strip(), m.group(6).strip())
            continue
        kind = m.group(1)
        args = [a.strip() for a in re.findall(r'"[^"]*"|[^,"]+', m.group(2)) if a.strip()]
        if kind == "AddMaterial":
            name = args[0].strip('"')
            if name == "G4_Galactic":
                dens = 1e-25                   # universe_mean_density in g/cm3
            else:
                dens = number(args[1])
            Z, pot = int(args[2]), number(args[3])
            ncomp = int(args[4]) if len(args) > 4 else 1
            state = STATES[args[5]] if len(args) > 5 else "solid"
            current = {"density": dens, "Z": Z, "I": pot, "state": state, "components": [],
                       "bycount": False, "formula": "", "gas": None}
            materials[name] = current
            if ncomp == 1 and Z > 0:
                current["components"].append((Z, 1.0))
                current["bycount"] = True
        else:
            el = args[0].strip('"')
            Z = zof[el] if not el.isdigit() else int(el)
            current["components"].append((Z, number(args[1])))
            if kind == "AddElementByAtomCount":
                current["bycount"] = True
    for mat in materials.values():
        if not mat["bycount"]:
            s = sum(w for _, w in mat["components"])
            mat["components"] = [(z, w / s) for z, w in mat["components"]]
    return materials


def main(elmfile, matfile):
    elements = parse_elements(open(elmfile).read())
    symbols = [(s, z) for s, z, _ in elements]
    materials = parse_materials(open(matfile).read(), symbols)
    simple = [n for n in materials if materials[n]["Z"] > 0 and n == "G4_" + symbols[materials[n]["Z"] - 1][0]]
    # the mean excitation energy of an element is the one of its simple material (G4NistManager)
    ion = {materials[n]["Z"]: materials[n]["I"] for n in simple}

    out = sys.stdout
    out.write("#  Generated by utils/generate_nist_data.py from the NIST builders of Geant4 %s\n" % GEANT4_VERSION)
    out.write("#  (source/materials/src/G4NistElementBuilder.cc and G4NistMaterialBuilder.cc). Do not edit.\n")
    out.write("#  Geant4 is distributed under the Geant4 Software License (https://geant4.web.cern.ch/license).\n\n")
    out.write("# Elements: (symbol, Z, mean excitation energy [eV], natural isotopes as (N, A [g/mole], abundance))\n")
    out.write("const _NIST_ELEMENTS = (\n")
    for sym, Z, isos in elements:
        isotxt = ", ".join("(%d, %r, %r)" % (n, a, w) for n, a, w in isos)
        out.write('    (symbol="%s", Z=%d, I=%r, isotopes=(%s%s)),\n' % (sym, Z, ion.get(Z, 0.0), isotxt, "," if len(isos) == 1 else ""))
    out.write(")\n\n")
    out.write("# Materials: density [g/cm3], mean excitation energy [eV] (0 if not defined), state, temperature\n")
    out.write("# [K] and pressure [atmosphere] when they are not the default ones, chemical formula, and the\n")
    out.write("# components as (Z, number of atoms) if bycount, else (Z, mass fraction)\n")
    out.write("const _NIST_MATERIALS = Dict{String,Any}(\n")
    names = simple + [n for n in MATERIALS if n not in simple]
    for n in names:
        if n not in materials:
            raise SystemExit("material %s not found in %s" % (n, matfile))
        m = materials[n]
        comps = ", ".join("(%d, %r)" % (z, w) for z, w in m["components"])
        if len(m["components"]) == 1:
            comps += ","
        T, P = "293.15", "1.0"
        if m["gas"] is not None:            # G4_Galactic: 2.73 kelvin, 3.e-18 hep_pascal
            T = number(m["gas"][0].split("*")[0])
            P = number(m["gas"][1].split("*")[0]) / 101325.0
        out.write('    "%s" => (density=%r, I=%r, state=:%s, temperature=%s, pressure=%s, formula="%s", bycount=%s, components=(%s)),\n'
                  % (n, m["density"], m["I"], m["state"], T, P, m["formula"], "true" if m["bycount"] else "false", comps))
    out.write(")\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
