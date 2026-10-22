"""Kessler warm-rain microphysics (moist idealized physics).

Kessler works in potential temperature and dry mixing ratios; CAM's state
is temperature and wet mixing ratios.  The parameterization therefore
converts on the way in and out around the Kessler scheme itself.
"""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization
from categories.state import geopotential_temp, qneg

calc_exner                     = ccpp_scheme_ref("calc_exner")
temp_to_potential_temp         = ccpp_scheme_ref("temp_to_potential_temp")
calc_dry_air_ideal_gas_density = ccpp_scheme_ref("calc_dry_air_ideal_gas_density")
wet_to_dry_water_vapor         = ccpp_scheme_ref("wet_to_dry_water_vapor")
wet_to_dry_cloud_liquid_water  = ccpp_scheme_ref("wet_to_dry_cloud_liquid_water")
wet_to_dry_rain                = ccpp_scheme_ref("wet_to_dry_rain")
kessler_scheme                 = ccpp_scheme_ref("kessler")
potential_temp_to_temp         = ccpp_scheme_ref("potential_temp_to_temp")
dry_to_wet_water_vapor         = ccpp_scheme_ref("dry_to_wet_water_vapor")
dry_to_wet_cloud_liquid_water  = ccpp_scheme_ref("dry_to_wet_cloud_liquid_water")
dry_to_wet_rain                = ccpp_scheme_ref("dry_to_wet_rain")
kessler_update                 = ccpp_scheme_ref("kessler_update")
kessler_diagnostics            = ccpp_scheme_ref("kessler_diagnostics")


@ccpp_composite
def kessler_microphysics():
    # In: temperature -> potential temperature, wet -> dry mixing ratios
    calc_exner()
    temp_to_potential_temp()
    calc_dry_air_ideal_gas_density()
    wet_to_dry_water_vapor()
    wet_to_dry_cloud_liquid_water()
    wet_to_dry_rain()

    kessler_scheme()

    # Out: back to temperature and wet mixing ratios, update the state
    potential_temp_to_temp()
    dry_to_wet_water_vapor()
    dry_to_wet_cloud_liquid_water()
    dry_to_wet_rain()
    kessler_update()
    qneg()
    geopotential_temp()


kessler = Parameterization(name="Kessler microphysics", run=kessler_microphysics,
                           diagnostics=kessler_diagnostics)
