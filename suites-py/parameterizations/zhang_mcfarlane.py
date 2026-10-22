"""Zhang-McFarlane deep convection (CAM's convect_deep).

ZM runs as four steps -- convection, evaporation of its precipitation,
momentum transport, constituent transport.  Each step computes tendencies,
writes their history output, and updates the state before the next step.

Two forms: ``zhang_mcfarlane`` (CAM4, CAM5) and
``zhang_mcfarlane_saving_tendency`` (CAM7), which also saves each step's
temperature tendency (save_ttend_from_convect_deep).
"""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization
from categories.state import (
    apply_constituent_tendencies,
    apply_heating_rate,
    apply_tendency_of_eastward_wind,
    apply_tendency_of_northward_wind,
    geopotential_temp,
    qneg,
)

zm_conv_options                  = ccpp_scheme_ref("zm_conv_options")
zm_convr                         = ccpp_scheme_ref("zm_convr")
zm_convr_tendency_diagnostics    = ccpp_scheme_ref("zm_convr_tendency_diagnostics")
zm_conv_evap                     = ccpp_scheme_ref("zm_conv_evap")
zm_evap_tendency_diagnostics     = ccpp_scheme_ref("zm_evap_tendency_diagnostics")
zm_conv_momtran                  = ccpp_scheme_ref("zm_conv_momtran")
zm_momtran_tendency_diagnostics  = ccpp_scheme_ref("zm_momtran_tendency_diagnostics")
zm_conv_convtran                 = ccpp_scheme_ref("zm_conv_convtran")
zm_tendency_diagnostics          = ccpp_scheme_ref("zm_tendency_diagnostics")
zm_diagnostics                   = ccpp_scheme_ref("zm_diagnostics")
zm_prepare_flux_for_check_energy = ccpp_scheme_ref("zm_prepare_flux_for_check_energy")
save_ttend_from_convect_deep     = ccpp_scheme_ref("save_ttend_from_convect_deep")
cloud_fraction_fice              = ccpp_scheme_ref("cloud_fraction_fice")
set_deep_conv_fluxes_to_general  = ccpp_scheme_ref("set_deep_conv_fluxes_to_general")
set_general_conv_fluxes_to_deep  = ccpp_scheme_ref("set_general_conv_fluxes_to_deep")


def _apply_heating_and_constituents():
    apply_heating_rate()
    apply_constituent_tendencies()
    qneg()
    geopotential_temp()


def _no_step():
    pass


def _zm_steps(save_tendency):
    """ZM's four steps; save_tendency runs after each of steps 1-3's
    diagnostics."""
    # 1. Convection
    zm_convr()
    zm_convr_tendency_diagnostics()
    save_tendency()
    _apply_heating_and_constituents()
    cloud_fraction_fice()

    # 2. Evaporation of convective precipitation (with the general
    #    convective-flux names zm_conv_evap reads)
    set_deep_conv_fluxes_to_general()
    zm_conv_evap()
    set_general_conv_fluxes_to_deep()
    zm_evap_tendency_diagnostics()
    save_tendency()
    _apply_heating_and_constituents()
    cloud_fraction_fice()

    # 3. Convective momentum transport
    zm_conv_momtran()
    zm_momtran_tendency_diagnostics()
    save_tendency()
    apply_heating_rate()
    apply_tendency_of_eastward_wind()
    apply_tendency_of_northward_wind()
    geopotential_temp()

    # 4. Convective constituent transport
    zm_conv_convtran()
    zm_tendency_diagnostics()
    apply_constituent_tendencies()
    qneg()
    geopotential_temp()

    zm_diagnostics()


@ccpp_composite
def zm_deep_convection():
    _zm_steps(_no_step)


@ccpp_composite
def zm_deep_convection_saving_tendency():
    _zm_steps(save_ttend_from_convect_deep)


zhang_mcfarlane = Parameterization(
    name="Zhang-McFarlane deep convection",
    setup=zm_conv_options,                                # ZM namelist (init only)
    run=zm_deep_convection,
    energy_flux_correction=zm_prepare_flux_for_check_energy)

zhang_mcfarlane_saving_tendency = Parameterization(
    name="Zhang-McFarlane deep convection, saving its temperature tendency",
    setup=zm_conv_options,
    run=zm_deep_convection_saving_tendency,
    energy_flux_correction=zm_prepare_flux_for_check_energy)
