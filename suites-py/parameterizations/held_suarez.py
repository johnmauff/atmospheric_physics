"""Held-Suarez (1994) idealized forcing: relaxation of temperature towards a
prescribed profile and Rayleigh friction of the winds (dry dynamical-core
test)."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization
from categories.state import (
    apply_heating_rate,
    apply_tendency_of_eastward_wind,
    apply_tendency_of_northward_wind,
    geopotential_temp,
)

held_suarez_1994_scheme = ccpp_scheme_ref("held_suarez_1994")


@ccpp_composite
def held_suarez_forcing():
    held_suarez_1994_scheme()
    apply_tendency_of_eastward_wind()
    apply_tendency_of_northward_wind()
    apply_heating_rate()
    geopotential_temp()


held_suarez = Parameterization(name="Held-Suarez (1994) forcing", run=held_suarez_forcing)
