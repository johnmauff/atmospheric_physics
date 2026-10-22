"""RRTMGP radiation as CAM-SIMA runs it: CAM's adapter around the
exchangeable RRTMGP package (rrtmgp_package.py, unchanged from xdsl-ccpp's
examples/rrtmgp_composite).

The adapter holds CAM policy that does not belong inside RRTMGP: when
shortwave and longwave run (``rrtmgp_pre``), where CAM's aerosol optics run,
which CAM schemes fill the package's hooks, and what happens after
radiation (heating carried between radiation steps, net heating, coupler
output, freeing RRTMGP's objects, the state update).
"""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization, nothing
from categories.state import update_temperature
from parameterizations.rrtmgp_package import rrtmgp_longwave, rrtmgp_prepare, rrtmgp_shortwave

rrtmgp_pre                        = ccpp_scheme_ref("rrtmgp_pre")
rrtmgp_dry_static_energy_tendency = ccpp_scheme_ref("rrtmgp_dry_static_energy_tendency")
calculate_net_heating             = ccpp_scheme_ref("calculate_net_heating")
rrtmgp_post                       = ccpp_scheme_ref("rrtmgp_post")
rrtmgp_diagnostics                = ccpp_scheme_ref("rrtmgp_diagnostics")

# CAM schemes that fill the package's hooks.
tropopause_find          = ccpp_scheme_ref("tropopause_find")
rrtmgp_constituents      = ccpp_scheme_ref("rrtmgp_constituents")
solar_irradiance_data    = ccpp_scheme_ref("solar_irradiance_data")
rrtmgp_cloud_diagnostics = ccpp_scheme_ref("rrtmgp_cloud_diagnostics")
rrtmgp_sw_diagnostics    = ccpp_scheme_ref("rrtmgp_sw_diagnostics")
rrtmgp_lw_diagnostics    = ccpp_scheme_ref("rrtmgp_lw_diagnostics")
rrtmgp_subcycle          = ccpp_scheme_ref("rrtmgp_subcycle")

# Radiation calls per step: one per gas set for forcing diagnostics, the
# last being the climate calculation.  A run-time count from CAM's namelist.
NDIAG = "number_of_diagnostic_subcycles"


@ccpp_composite
def cam_rrtmgp(*, aerosol_optics=nothing):
    """aerosol_optics runs right after rrtmgp_pre, which provides the night
    columns its diagnostics need."""
    rrtmgp_pre()                                  # when SW/LW run this step
    aerosol_optics()
    rrtmgp_prepare(tropopause=tropopause_find)
    rrtmgp_shortwave(n_calls=NDIAG,
                     next_call=rrtmgp_subcycle,
                     gas_source=rrtmgp_constituents,
                     solar_source=solar_irradiance_data,
                     aerosol_optics=nothing,      # CAM's run before the package
                     cloud_diagnostics=rrtmgp_cloud_diagnostics,
                     diagnostics=rrtmgp_sw_diagnostics)
    rrtmgp_longwave(n_calls=NDIAG,
                    next_call=rrtmgp_subcycle,
                    gas_source=rrtmgp_constituents,
                    diagnostics=rrtmgp_lw_diagnostics)
    rrtmgp_dry_static_energy_tendency()           # heating carried between radiation steps
    calculate_net_heating()                       # net radiative heating and fluxes
    rrtmgp_post()                                 # coupler output, free RRTMGP objects
    rrtmgp_diagnostics()                          # history output
    update_temperature()


rrtmgp = Parameterization(name="RRTMGP", run=cam_rrtmgp)
