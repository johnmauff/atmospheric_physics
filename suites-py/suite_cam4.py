"""CAM4 physics suite: Python version of suites/suite_cam4.xml.

A suite is a configuration: each line is a physics *category* (from
categories/), and its arguments are the *parameterizations* (from
parameterizations/) that fill the category's roles in CAM4.  What a
category does with its roles -- order, energy checks, state updates -- is
defined once in the category and shared by every CAM version.
"""

from xdsl_ccpp.frontend.py_api import ccpp_suite, export_suite

from categories.clouds import clouds
from categories.convection import convection
from categories.diagnostics import (
    state_diagnostics,
    surface_coupling,
    tendency_diagnostics,
    tropopause_output,
)
from categories.energy import end_of_step_energy, energy_fixer
from categories.gravity_wave_drag import gravity_wave_drag
from categories.inputs import prescribed_inputs
from categories.radiation import radiation
from categories.vertical_diffusion import vertical_diffusion
from parameterizations.aerosol_optics import bulk_aerosol_optics
from parameterizations.dry_adiabatic_adjustment import dadadj
from parameterizations.gravity_wave_sources import orographic
from parameterizations.hack import hack
from parameterizations.holtslag_boville import holtslag_boville
from parameterizations.prescribed import (
    prescribed_aerosols_and_deposition,
    prescribed_ozone,
    prescribed_volcanic_aerosols,
)
from parameterizations.rasch_kristjansson import rasch_kristjansson
from parameterizations.rrtmgp import rrtmgp
from parameterizations.zhang_mcfarlane import zhang_mcfarlane


# In a CCPP suite class, each method is a group recorded by @ccpp_suite,
# not a real method: it takes no 'self', and the class holds only groups.
# pylint: disable=no-method-argument,too-few-public-methods
@ccpp_suite("cam4", version="1.0")
class Cam4Suite:
    """The 'cam4' suite."""

    def physics_before_coupler():
        """Physics before the coupler."""
        prescribed_inputs(ozone=prescribed_ozone,
                          aerosols=prescribed_aerosols_and_deposition,
                          volcanic_aerosols=prescribed_volcanic_aerosols)
        energy_fixer()
        convection(dry_adjustment=dadadj, deep=zhang_mcfarlane, shallow=hack)
        clouds(macro_micro=rasch_kristjansson)
        state_diagnostics()
        radiation(scheme=rrtmgp, aerosol_optics=bulk_aerosol_optics)
        tropopause_output()
        surface_coupling()

    def physics_after_coupler():
        """Physics after the coupler.  (No chemistry or aerosol processes.)"""
        vertical_diffusion(boundary_layer=holtslag_boville)
        gravity_wave_drag(sources=[orographic])
        # Not yet: thermo_water_update, dme_adjust (atmospheric_physics#222).
        end_of_step_energy()
        tendency_diagnostics()


if __name__ == "__main__":
    export_suite(Cam4Suite)
