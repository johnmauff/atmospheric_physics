"""Category: history output that belongs to no single process, and the
surface variables for the coupler."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

sima_state_diagnostics         = ccpp_scheme_ref("sima_state_diagnostics")
sima_tend_diagnostics          = ccpp_scheme_ref("sima_tend_diagnostics")
tropopause_find                = ccpp_scheme_ref("tropopause_find")
tropopause_diagnostics         = ccpp_scheme_ref("tropopause_diagnostics")
calc_dry_air_ideal_gas_density = ccpp_scheme_ref("calc_dry_air_ideal_gas_density")
set_surface_coupling_vars      = ccpp_scheme_ref("set_surface_coupling_vars")


@ccpp_composite
def state_diagnostics(*, parameterizations=()):
    """History output of the physics state, then of *parameterizations*
    (each one's ``diagnostics``) that write theirs alongside it."""
    sima_state_diagnostics()
    for parameterization in parameterizations:
        parameterization.diagnostics()


@ccpp_composite
def tendency_diagnostics():
    """History output of the physics tendencies (end of the step)."""
    sima_tend_diagnostics()


@ccpp_composite
def tropopause_output():
    """Find the tropopause and write it out.  Runs after convective
    adjustment, so supersaturated cells are not written to the tape."""
    tropopause_find()
    tropopause_diagnostics()


@ccpp_composite
def surface_coupling():
    """Surface variables the coupler reads.  Must be the last step before
    the coupler."""
    calc_dry_air_ideal_gas_density()
    set_surface_coupling_vars()
