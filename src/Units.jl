"""
    Geom4hep.Units

System of units, with the conventions of Geant4 (CLHEP) for lengths (mm), angles (radian), energies
(MeV) and temperatures (kelvin). Masses are in mg, so that densities are in g/cm³ (`g/cm3 == 1`)
and molar masses in mg/mole (`g/mole == 1000`). Use as `using Geom4hep.Units;
Box{Float64}(1cm, 2cm, 3cm)` or `density = 7.87g/cm3`.

It also defines the physical constants used for the materials: `Avogadro`, `k_Boltzmann`, `amu`,
`fine_structure_const`, `classic_electr_radius`, the standard temperature and pressure
(`STP_Temperature`, `NTP_Temperature`, `STP_Pressure`) and `universe_mean_density`, with the
values of CLHEP.
"""
module Units

export mm, cm, m, km, um, nm, mm2, cm2, m2, mm3, cm3, m3, radian, rad, mrad, degree, deg,
       mg, g, kg, mole, eV, keV, MeV, GeV, TeV, kelvin, pascal, hep_pascal, bar, atmosphere

#--- Length (internal unit: mm)
const mm = 1.0
const cm = 10mm
const m  = 1000mm
const km = 1000m
const um = 1e-3mm
const nm = 1e-6mm
const mm2 = mm^2
const cm2 = cm^2
const m2  = m^2
const mm3 = mm^3
const cm3 = cm^3
const m3  = m^3

#--- Angle (internal unit: radian)
const radian = 1.0
const rad = radian
const mrad = 1e-3rad
const degree = π / 180.0 * radian
const deg = degree

#--- Mass (internal unit: mg, so that g/cm3 = 1) and amount of substance
const mg = 1.0
const g  = 1000mg
const kg = 1000g
const mole = 1.0

#--- Energy (internal unit: MeV)
const MeV = 1.0
const eV  = 1e-6MeV
const keV = 1e-3MeV
const GeV = 1e3MeV
const TeV = 1e6MeV

#--- Temperature (internal unit: kelvin) and pressure (MeV/mm³ as in CLHEP)
const kelvin = 1.0
const _e_SI = 1.602176634e-19                     # positron charge in coulomb (= joule/eV)
const pascal = 1e-6MeV / _e_SI / m3              # joule/m3
const hep_pascal = pascal
const bar = 1e5pascal
const atmosphere = 101325pascal

#--- Physical constants (CLHEP/Units/PhysicalConstants.h of Geant4 11)
const Avogadro = 6.02214076e+23 / mole
const k_Boltzmann = 8.617333e-11MeV / kelvin
const STP_Temperature = 273.15kelvin
const NTP_Temperature = 293.15kelvin
const STP_Pressure = 1atmosphere
const universe_mean_density = 1e-25g / cm3

const _c_light = 2.99792458e+8                      # m/s
const _h_Planck = 6.62607015e-34                    # J s
const _mu0 = 4π * 1e-7                              # henry/m
const electron_mass_c2 = 0.510998910MeV
const amu_c2 = 931.494028MeV
"Atomic mass unit (as a mass, CLHEP `amu = amu_c2/c_squared`)."
const amu = amu_c2 / MeV * 1e6 * _e_SI / _c_light^2 * kg
"Fine structure constant (CLHEP `elm_coupling/hbarc`)."
const fine_structure_const = _e_SI^2 * _c_light * _mu0 / (2 * _h_Planck)
const hbarc = _h_Planck * _c_light / 2π / _e_SI * eV * m
"Classical electron radius (CLHEP `elm_coupling/electron_mass_c2`)."
const classic_electr_radius = fine_structure_const * hbarc / electron_mass_c2

end # module Units
