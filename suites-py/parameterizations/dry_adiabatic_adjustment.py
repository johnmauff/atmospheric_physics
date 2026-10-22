"""Dry adiabatic adjustment (dadadj): removes dry static instability."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization
from categories.state import apply_constituent_tendencies, apply_heating_rate, geopotential_temp, qneg

dadadj_scheme = ccpp_scheme_ref("dadadj")


@ccpp_composite
def dry_adiabatic_adjustment():
    dadadj_scheme()
    apply_constituent_tendencies()
    apply_heating_rate()
    qneg()
    geopotential_temp()


dadadj = Parameterization(name="dry adiabatic adjustment", run=dry_adiabatic_adjustment)
