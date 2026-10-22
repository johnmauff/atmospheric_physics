"""Beljaars orographic form drag from subgrid orography (CAM5, CAM7).

Contributes to two categories: the form drag itself (surface drag), and an
addition to the wind damping rate inside vertical diffusion.
"""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization

beljaars_drag                  = ccpp_scheme_ref("beljaars_drag")
beljaars_drag_diagnostics      = ccpp_scheme_ref("beljaars_drag_diagnostics")
beljaars_add_wind_damping_rate = ccpp_scheme_ref("beljaars_add_wind_damping_rate")


@ccpp_composite
def beljaars_form_drag():
    beljaars_drag()
    beljaars_drag_diagnostics()


beljaars = Parameterization(name="Beljaars form drag", run=beljaars_form_drag)
beljaars_wind_damping = Parameterization(name="Beljaars wind damping",
                                         run=beljaars_add_wind_damping_rate)
