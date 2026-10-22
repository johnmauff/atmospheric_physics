"""Thatcher-Jablonowski (2016) moist idealized physics: large-scale
precipitation, and a surface flux / boundary layer scheme with Held-Suarez
temperature relaxation."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization
from categories.state import (
    apply_heating_rate,
    apply_tendency_of_eastward_wind,
    apply_tendency_of_northward_wind,
    qneg,
)

tj2016_precip_scheme     = ccpp_scheme_ref("tj2016_precip")
tj2016_sfc_pbl_hs_scheme = ccpp_scheme_ref("tj2016_sfc_pbl_hs")


@ccpp_composite
def tj2016_large_scale_precipitation():
    tj2016_precip_scheme()
    apply_heating_rate()
    qneg()


@ccpp_composite
def tj2016_surface_and_boundary_layer():
    tj2016_sfc_pbl_hs_scheme()
    apply_heating_rate()
    apply_tendency_of_eastward_wind()
    apply_tendency_of_northward_wind()
    qneg()


tj2016_precipitation = Parameterization(name="TJ2016 large-scale precipitation",
                                        run=tj2016_large_scale_precipitation)
tj2016_surface_pbl = Parameterization(name="TJ2016 surface fluxes and boundary layer",
                                      run=tj2016_surface_and_boundary_layer)
