"""Hack shallow convection (CAM4).

Runs the shallow scheme, then evaporates its precipitation with ZM's
evaporation scheme, then sums deep and shallow quantities (total convective
rainout).  CAM's convect_shallow takes in total precipitation from both, so
all three steps are one process for the energy check.
"""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization
from categories.state import apply_constituent_tendencies, apply_heating_rate, geopotential_temp, qneg

convect_shallow_diagnostics = ccpp_scheme_ref("convect_shallow_diagnostics")
hack_convect_shallow        = ccpp_scheme_ref("hack_convect_shallow")
convect_shallow_diagnostics_after_shallow_scheme = ccpp_scheme_ref(
    "convect_shallow_diagnostics_after_shallow_scheme")
cloud_fraction_fice                = ccpp_scheme_ref("cloud_fraction_fice")
set_shallow_conv_fluxes_to_general = ccpp_scheme_ref("set_shallow_conv_fluxes_to_general")
zm_conv_evap                       = ccpp_scheme_ref("zm_conv_evap")
set_general_conv_fluxes_to_shallow = ccpp_scheme_ref("set_general_conv_fluxes_to_shallow")
convect_shallow_diagnostics_after_convective_evaporation = ccpp_scheme_ref(
    "convect_shallow_diagnostics_after_convective_evaporation")
convect_shallow_sum_to_deep = ccpp_scheme_ref("convect_shallow_sum_to_deep")
convect_shallow_diagnostics_after_sum_to_deep = ccpp_scheme_ref(
    "convect_shallow_diagnostics_after_sum_to_deep")


def _apply_heating_and_constituents():
    apply_heating_rate()
    apply_constituent_tendencies()
    qneg()
    geopotential_temp()


@ccpp_composite
def hack_shallow_convection():
    hack_convect_shallow()
    convect_shallow_diagnostics_after_shallow_scheme()
    _apply_heating_and_constituents()
    cloud_fraction_fice()

    # Subcloud evaporation: shallow outputs renamed to the generic inputs
    # zm_conv_evap reads, and back.
    set_shallow_conv_fluxes_to_general()
    zm_conv_evap()
    set_general_conv_fluxes_to_shallow()
    convect_shallow_diagnostics_after_convective_evaporation()
    _apply_heating_and_constituents()

    # Sum deep convective quantities with shallow (total rainout).
    convect_shallow_sum_to_deep()
    convect_shallow_diagnostics_after_sum_to_deep()


hack = Parameterization(
    name="Hack shallow convection",
    setup=convect_shallow_diagnostics,      # registers its history fields (init only)
    run=hack_shallow_convection)
