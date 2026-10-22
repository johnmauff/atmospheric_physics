"""State-update schemes: applying a step's tendencies to the physics state.

After a step computes tendencies, CAM applies them -- temperature (through
the heating rate), constituents (then clipping negative values), winds --
and recomputes geopotential.  The order differs between parameterizations
(heating before or after constituents; north or east wind first), so each
parameterization writes its own update with these schemes.  Only sequences
that recur unchanged get a helper.
"""

from xdsl_ccpp.frontend.py_api import ccpp_scheme_ref

apply_heating_rate                = ccpp_scheme_ref("apply_heating_rate")
apply_constituent_tendencies      = ccpp_scheme_ref("apply_constituent_tendencies")
apply_tendency_of_eastward_wind   = ccpp_scheme_ref("apply_tendency_of_eastward_wind")
apply_tendency_of_northward_wind  = ccpp_scheme_ref("apply_tendency_of_northward_wind")
apply_tendency_of_air_temperature = ccpp_scheme_ref("apply_tendency_of_air_temperature")
qneg                              = ccpp_scheme_ref("qneg")
geopotential_temp                 = ccpp_scheme_ref("geopotential_temp")
update_dry_static_energy          = ccpp_scheme_ref("update_dry_static_energy")


def update_temperature():
    """Apply the heating rate, then recompute geopotential."""
    apply_heating_rate()
    geopotential_temp()
