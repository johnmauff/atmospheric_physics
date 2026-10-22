"""CAM5 physics suite: Python version of suites/suite_cam5.xml.

A suite is a configuration: each line is a physics *category* (from
categories/), and its arguments are the *parameterizations* (from
parameterizations/) that fill the category's roles in CAM5.

Not implemented yet: UW shallow convection, Park macrophysics (a cloud
fraction placeholder runs instead), chemistry and aerosol processes.
"""

from dataclasses import replace

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
from categories.surface_drag import surface_drag
from categories.vertical_diffusion import vertical_diffusion
from parameterizations.beljaars import beljaars, beljaars_wind_damping
from parameterizations.bretherton_park import uw_bretherton_park
from parameterizations.cloud_water import cloud_water_paths
from parameterizations.dropmixnuc import dropmixnuc_surface_fluxes
from parameterizations.dry_adiabatic_adjustment import dadadj
from parameterizations.gravity_wave_sources import orographic_without_diagnostics
from parameterizations.microp_aero import microp_aero_tke
from parameterizations.park_macrophysics import park_placeholder
from parameterizations.prescribed import prescribed_bulk_aerosols
from parameterizations.pumas import pumas
from parameterizations.rrtmgp import rrtmgp
from parameterizations.turbulent_mountain_stress import turbulent_mountain_stress
from parameterizations.zhang_mcfarlane import zhang_mcfarlane

# CAM5's suite checks ZM's energy budget without ZM's flux correction (CAM4
# and CAM7 apply it; likely an omission in suite_cam5.xml).
zhang_mcfarlane_cam5 = replace(zhang_mcfarlane, energy_flux_correction=None)


# In a CCPP suite class, each method is a group recorded by @ccpp_suite,
# not a real method: it takes no 'self', and the class holds only groups.
# pylint: disable=no-method-argument,too-few-public-methods
@ccpp_suite("cam5", version="1.0")
class Cam5Suite:
    """The 'cam5' suite."""

    def physics_before_coupler():
        """Physics before the coupler."""
        prescribed_inputs(aerosols=prescribed_bulk_aerosols)     # for the microphysics
        energy_fixer()
        convection(dry_adjustment=dadadj, deep=zhang_mcfarlane_cam5)
        clouds(macrophysics=park_placeholder,
               aerosol_activation=microp_aero_tke,
               microphysics=pumas)
        state_diagnostics()
        radiation(scheme=rrtmgp, cloud_water=cloud_water_paths)
        tropopause_output()
        surface_coupling()

    def physics_after_coupler():
        """Physics after the coupler."""
        surface_drag(form_drag=beljaars)
        vertical_diffusion(boundary_layer=uw_bretherton_park,
                           mountain_stress=turbulent_mountain_stress,
                           form_drag_damping=beljaars_wind_damping,
                           undiffused_tracer_fluxes=dropmixnuc_surface_fluxes)
        gravity_wave_drag(sources=[orographic_without_diagnostics])
        # Not yet: dme_adjust (atmospheric_physics#222).
        end_of_step_energy()
        tendency_diagnostics()


if __name__ == "__main__":
    export_suite(Cam5Suite)
