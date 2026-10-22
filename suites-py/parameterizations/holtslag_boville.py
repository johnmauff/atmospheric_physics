"""Holtslag-Boville boundary layer (CAM4): diffusion coefficients for
vertical diffusion."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization

holtslag_boville_diff_options       = ccpp_scheme_ref("holtslag_boville_diff_options")
holtslag_boville_diff               = ccpp_scheme_ref("holtslag_boville_diff")
hb_pbl_independent_coefficients     = ccpp_scheme_ref("hb_pbl_independent_coefficients")
compute_kinematic_fluxes_and_obklen = ccpp_scheme_ref("compute_kinematic_fluxes_and_obklen")
hb_pbl_dependent_coefficients       = ccpp_scheme_ref("hb_pbl_dependent_coefficients")
hb_diff_exchange_coefficients       = ccpp_scheme_ref("hb_diff_exchange_coefficients")
holtslag_boville_diff_diagnostics   = ccpp_scheme_ref("holtslag_boville_diff_diagnostics")


@ccpp_composite
def holtslag_boville_coefficients():
    holtslag_boville_diff()                 # init only
    hb_pbl_independent_coefficients()       # boundary-layer-independent coefficients
    compute_kinematic_fluxes_and_obklen()
    hb_pbl_dependent_coefficients()         # boundary-layer-dependent coefficients
    hb_diff_exchange_coefficients()


holtslag_boville = Parameterization(
    name="Holtslag-Boville boundary layer",
    setup=holtslag_boville_diff_options,    # namelist
    run=holtslag_boville_coefficients,
    diagnostics=holtslag_boville_diff_diagnostics)
