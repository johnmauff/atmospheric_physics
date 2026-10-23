"""Kessler warm-rain microphysics (moist idealized physics) as CAM runs it:
CAM's adapter around the Kessler package (kessler_package.py).

Kessler works in dry mixing ratios; CAM's are wet. The adapter fills the
package's conversion hooks with CAM's schemes for the three water species,
then cleans up negative mixing ratios and recomputes the geopotential.
"""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization
from categories.state import geopotential_temp, qneg
from parameterizations.kessler_package import kessler_microphysics

wet_to_dry_water_vapor         = ccpp_scheme_ref("wet_to_dry_water_vapor")
wet_to_dry_cloud_liquid_water  = ccpp_scheme_ref("wet_to_dry_cloud_liquid_water")
wet_to_dry_rain                = ccpp_scheme_ref("wet_to_dry_rain")
dry_to_wet_water_vapor         = ccpp_scheme_ref("dry_to_wet_water_vapor")
dry_to_wet_cloud_liquid_water  = ccpp_scheme_ref("dry_to_wet_cloud_liquid_water")
dry_to_wet_rain                = ccpp_scheme_ref("dry_to_wet_rain")
kessler_diagnostics            = ccpp_scheme_ref("kessler_diagnostics")


@ccpp_composite
def water_to_dry():
    """Kessler's water species, wet to dry mixing ratios."""
    wet_to_dry_water_vapor()
    wet_to_dry_cloud_liquid_water()
    wet_to_dry_rain()


@ccpp_composite
def water_to_wet():
    """Kessler's water species, dry to wet mixing ratios."""
    dry_to_wet_water_vapor()
    dry_to_wet_cloud_liquid_water()
    dry_to_wet_rain()


@ccpp_composite
def cam_kessler():
    kessler_microphysics(to_dry=water_to_dry, to_wet=water_to_wet)
    qneg()
    geopotential_temp()


kessler = Parameterization(name="Kessler microphysics", run=cam_kessler,
                           diagnostics=kessler_diagnostics)
