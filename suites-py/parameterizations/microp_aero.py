"""microp_aero: coupling between aerosols and cloud microphysics (from CAM5),
with nochem BAM (bulk) aerosols.

Computes the subgrid vertical velocity, then ice nucleation and droplet
activation, which the microphysics (PUMAS) uses.  Model versions differ in
where the subgrid vertical velocity comes from: TKE in CAM5, CLUBB in CAM7.
"""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization

tropopause_find                             = ccpp_scheme_ref("tropopause_find")
compute_subgrid_vertical_velocity_tke       = ccpp_scheme_ref("compute_subgrid_vertical_velocity_tke")
compute_subgrid_vertical_velocity_clubb     = ccpp_scheme_ref("compute_subgrid_vertical_velocity_clubb")
scale_subgrid_vertical_velocity             = ccpp_scheme_ref("scale_subgrid_vertical_velocity")
scale_subgrid_vertical_velocity_diagnostics = ccpp_scheme_ref("scale_subgrid_vertical_velocity_diagnostics")
dust_default_radii                          = ccpp_scheme_ref("dust_default_radii")
nucleate_ice_ccpp                           = ccpp_scheme_ref("nucleate_ice_ccpp")
nucleate_ice_diagnostics                    = ccpp_scheme_ref("nucleate_ice_diagnostics")
ndrop_bam_ccpp                              = ccpp_scheme_ref("ndrop_bam_ccpp")
ndrop_bam_diagnostics                       = ccpp_scheme_ref("ndrop_bam_diagnostics")
hetfrz_classnuc_stub                        = ccpp_scheme_ref("hetfrz_classnuc_stub")


def _microp_aero(subgrid_vertical_velocity):
    tropopause_find()                       # chemical tropopause
    subgrid_vertical_velocity()
    scale_subgrid_vertical_velocity()
    scale_subgrid_vertical_velocity_diagnostics()
    dust_default_radii()
    nucleate_ice_ccpp()                     # ice nucleation (Meyers 1992 without aerosol models)
    nucleate_ice_diagnostics()
    ndrop_bam_ccpp()                        # droplet activation
    ndrop_bam_diagnostics()
    hetfrz_classnuc_stub()                  # passes the hetfrz namelist option on


@ccpp_composite
def microp_aero_tke_velocity():
    """Subgrid vertical velocity from TKE."""
    _microp_aero(compute_subgrid_vertical_velocity_tke)


@ccpp_composite
def microp_aero_clubb_velocity():
    """Subgrid vertical velocity from CLUBB's TKE (WP2_nadv)."""
    _microp_aero(compute_subgrid_vertical_velocity_clubb)


microp_aero_tke = Parameterization(name="microp_aero (TKE subgrid velocity)",
                                   run=microp_aero_tke_velocity)
microp_aero_clubb = Parameterization(name="microp_aero (CLUBB subgrid velocity)",
                                     run=microp_aero_clubb_velocity)
