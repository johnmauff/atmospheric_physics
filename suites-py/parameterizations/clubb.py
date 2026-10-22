"""CLUBB (CAM7 macrophysics and shallow convection) -- not implemented yet.

``clubb_stub`` only provides the do_clubb_sgs flag the PUMAS interface
reads.  In CAM7, CLUBB, microp_aero and PUMAS run together inside the
macro/microphysics substep loop (cld_macmic_num_steps), which stays off
until CLUBB is implemented.
"""

from xdsl_ccpp.frontend.py_api import ccpp_scheme_ref

from categories.interface import Parameterization

clubb_stub                        = ccpp_scheme_ref("clubb_stub")
set_cloud_fraction_top            = ccpp_scheme_ref("set_cloud_fraction_top")
compute_cloud_fraction_two_moment = ccpp_scheme_ref("compute_cloud_fraction_two_moment")


def _cloud_fraction_setup():
    """Cloud fraction settings for CAM7's two-moment clouds (init only)."""
    set_cloud_fraction_top()                # top level for cloud fraction
    compute_cloud_fraction_two_moment()


clubb_placeholder = Parameterization(name="CLUBB placeholder (clubb_stub)",
                                     setup=_cloud_fraction_setup,
                                     run=clubb_stub)
