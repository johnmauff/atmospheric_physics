"""Droplet mixing and nucleation (dropmixnuc) -- here only its surface-flux
step for tracers vertical diffusion does not diffuse (CAM5)."""

from xdsl_ccpp.frontend.py_api import ccpp_scheme_ref

from categories.interface import Parameterization

dropmixnuc_apply_surface_fluxes = ccpp_scheme_ref("dropmixnuc_apply_surface_fluxes")

dropmixnuc_surface_fluxes = Parameterization(
    name="dropmixnuc surface fluxes for undiffused tracers", run=dropmixnuc_apply_surface_fluxes)
