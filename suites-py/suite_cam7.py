"""CAM7 physics suite (-chem none variant): Python version of
suites/suite_cam7.xml.

A suite is a configuration: each line is a physics *category* (from
categories/), and its arguments are the *parameterizations* (from
parameterizations/) that fill the category's roles in CAM7.

CAM7 moves clouds and radiation after the coupler.  Not implemented yet:
CLUBB (a stub stands in; it would also be CAM7's shallow convection and
vertical diffusion), the macro/microphysics substep loop, chemistry and
aerosol processes, the surface coupling variables.
"""

from xdsl_ccpp.frontend.py_api import ccpp_suite, export_suite

from categories.clouds import clouds
from categories.convection import convection
from categories.diagnostics import state_diagnostics, tendency_diagnostics, tropopause_output
from categories.energy import end_of_step_energy, energy_fixer
from categories.gravity_wave_drag import gravity_wave_drag
from categories.inputs import prescribed_inputs
from categories.radiation import radiation
from categories.surface_drag import surface_drag
from parameterizations.aerosol_optics import bulk_aerosol_optics_with_humidity
from parameterizations.beljaars import beljaars
from parameterizations.cloud_water import cloud_water_paths
from parameterizations.clubb import clubb_placeholder
from parameterizations.dry_adiabatic_adjustment import dadadj
from parameterizations.gravity_wave_sources import deep_convection, frontogenesis, moving_mountain, ridge
from parameterizations.microp_aero import microp_aero_clubb
from parameterizations.prescribed import (
    prescribed_aerosols_and_deposition,
    prescribed_ozone,
    prescribed_volcanic_aerosols,
)
from parameterizations.pumas import pumas
from parameterizations.rrtmgp import rrtmgp
from parameterizations.zhang_mcfarlane import zhang_mcfarlane_saving_tendency


# In a CCPP suite class, each method is a group recorded by @ccpp_suite,
# not a real method: it takes no 'self', and the class holds only groups.
# pylint: disable=no-method-argument,too-few-public-methods
@ccpp_suite("cam7", version="1.0")
class Cam7Suite:
    """The 'cam7' suite."""

    def physics_before_coupler():
        """Physics before the coupler."""
        prescribed_inputs(ozone=prescribed_ozone,
                          aerosols=prescribed_aerosols_and_deposition,
                          volcanic_aerosols=prescribed_volcanic_aerosols)
        energy_fixer(diagnostics=True)
        convection(dry_adjustment=dadadj, deep=zhang_mcfarlane_saving_tendency)
        # Not yet: surface_coupling(), once all other CAM7 schemes are present.

    def physics_after_coupler():
        """Physics after the coupler."""
        # In CAM7 this repeats cld_macmic_num_steps times (off until CLUBB exists).
        clouds(macrophysics=clubb_placeholder,
               aerosol_activation=microp_aero_clubb,
               microphysics=pumas)
        state_diagnostics()
        radiation(scheme=rrtmgp,
                  cloud_water=cloud_water_paths,
                  aerosol_optics=bulk_aerosol_optics_with_humidity)
        tropopause_output()
        surface_drag(form_drag=beljaars)
        # Vertical diffusion: CLUBB (not implemented).
        gravity_wave_drag(sources=[moving_mountain, deep_convection, frontogenesis, ridge])
        # Not yet: thermo_water_update, dme_adjust (atmospheric_physics#222).
        end_of_step_energy()
        tendency_diagnostics()


if __name__ == "__main__":
    export_suite(Cam7Suite)
