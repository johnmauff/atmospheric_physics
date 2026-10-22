"""UW (Bretherton-Park) eddy diffusivity boundary layer (CAM5; formerly
diag_TKE): diffusion coefficients for vertical diffusion."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization

bretherton_park_diff                  = ccpp_scheme_ref("bretherton_park_diff")
bretherton_park_diff_diagnostics      = ccpp_scheme_ref("bretherton_park_diff_diagnostics")
compute_kinematic_fluxes_and_obklen   = ccpp_scheme_ref("compute_kinematic_fluxes_and_obklen")
eddy_diffusivity_adjustment_above_pbl = ccpp_scheme_ref("eddy_diffusivity_adjustment_above_pbl")


@ccpp_composite
def bretherton_park_coefficients():
    bretherton_park_diff()
    bretherton_park_diff_diagnostics()
    compute_kinematic_fluxes_and_obklen()
    eddy_diffusivity_adjustment_above_pbl()


uw_bretherton_park = Parameterization(name="UW (Bretherton-Park) boundary layer",
                                      run=bretherton_park_coefficients)
