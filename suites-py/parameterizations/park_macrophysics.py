"""Park macrophysics (CAM5) -- not implemented yet.

CAM5's macrop_driver/cldwat2m is not CCPP-ized; the suite runs only the
cloud fraction calculation in its place.
"""

from xdsl_ccpp.frontend.py_api import ccpp_scheme_ref

from categories.interface import Parameterization

compute_cloud_fraction_two_moment = ccpp_scheme_ref("compute_cloud_fraction_two_moment")
compute_cloud_fraction            = ccpp_scheme_ref("compute_cloud_fraction")

park_placeholder = Parameterization(
    name="Park macrophysics placeholder (cloud fraction only)",
    setup=compute_cloud_fraction_two_moment,   # cloud fraction for two-moment microphysics (init only)
    run=compute_cloud_fraction)
