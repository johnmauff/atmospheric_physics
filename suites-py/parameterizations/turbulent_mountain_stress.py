"""Turbulent mountain stress (CAM5): adds to the surface drag coefficient in
vertical diffusion."""

from xdsl_ccpp.frontend.py_api import ccpp_scheme_ref

from categories.interface import Parameterization

turbulent_mountain_stress_add_drag_coefficient = ccpp_scheme_ref(
    "turbulent_mountain_stress_add_drag_coefficient")

turbulent_mountain_stress = Parameterization(
    name="turbulent mountain stress", run=turbulent_mountain_stress_add_drag_coefficient)
