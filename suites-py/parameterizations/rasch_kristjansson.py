"""Rasch-Kristjansson stratiform macro- and microphysics (CAM4).

Written as CAM's steps RK1-RK6.  Steps RK1-RK3 (sedimentation through
prognostic cloud water) are energy-checked as one block; the reserved
convective liquid detrained into the column is added to the boundary fluxes
before the check.
"""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.energy import energy_checked
from categories.interface import Parameterization
from categories.state import apply_constituent_tendencies, apply_heating_rate, geopotential_temp, qneg

tropopause_find                    = ccpp_scheme_ref("tropopause_find")
rk_stratiform_diagnostics          = ccpp_scheme_ref("rk_stratiform_diagnostics")
rk_stratiform_check_qtlcwat        = ccpp_scheme_ref("rk_stratiform_check_qtlcwat")
cloud_particle_sedimentation       = ccpp_scheme_ref("cloud_particle_sedimentation")
cloud_particle_sedimentation_diagnostics = ccpp_scheme_ref("cloud_particle_sedimentation_diagnostics")
rk_stratiform_sedimentation        = ccpp_scheme_ref("rk_stratiform_sedimentation")
rk_stratiform_detrain_convective_condensate = ccpp_scheme_ref("rk_stratiform_detrain_convective_condensate")
convective_cloud_cover             = ccpp_scheme_ref("convective_cloud_cover")
convective_cloud_cover_diagnostics = ccpp_scheme_ref("convective_cloud_cover_diagnostics")
compute_cloud_fraction             = ccpp_scheme_ref("compute_cloud_fraction")
compute_cloud_fraction_diagnostics = ccpp_scheme_ref("compute_cloud_fraction_diagnostics")
rk_stratiform_cloud_fraction_perturbation = ccpp_scheme_ref("rk_stratiform_cloud_fraction_perturbation")
rk_stratiform_cloud_fraction_perturbation_diagnostics = ccpp_scheme_ref(
    "rk_stratiform_cloud_fraction_perturbation_diagnostics")
rk_stratiform_external_forcings    = ccpp_scheme_ref("rk_stratiform_external_forcings")
cloud_fraction_fice                = ccpp_scheme_ref("cloud_fraction_fice")
prognostic_cloud_water             = ccpp_scheme_ref("prognostic_cloud_water")
rk_stratiform_condensate_repartioning = ccpp_scheme_ref("rk_stratiform_condensate_repartioning")
rk_stratiform_condensate_repartioning_diagnostics = ccpp_scheme_ref(
    "rk_stratiform_condensate_repartioning_diagnostics")
rk_stratiform_prognostic_cloud_water_tendencies = ccpp_scheme_ref(
    "rk_stratiform_prognostic_cloud_water_tendencies")
rk_stratiform_prognostic_cloud_water_tendencies_diagnostics = ccpp_scheme_ref(
    "rk_stratiform_prognostic_cloud_water_tendencies_diagnostics")
rk_stratiform_prepare_flux_for_check_energy = ccpp_scheme_ref("rk_stratiform_prepare_flux_for_check_energy")
rk_stratiform_cloud_optical_properties = ccpp_scheme_ref("rk_stratiform_cloud_optical_properties")
rk_stratiform_cloud_optical_properties_diagnostics = ccpp_scheme_ref(
    "rk_stratiform_cloud_optical_properties_diagnostics")
rk_stratiform_save_qtlcwat         = ccpp_scheme_ref("rk_stratiform_save_qtlcwat")


def _rk1_to_rk3():
    # RK1: cloud sedimentation (pcwsediment)
    cloud_particle_sedimentation()
    cloud_particle_sedimentation_diagnostics()
    apply_constituent_tendencies()
    apply_heating_rate()
    qneg()
    geopotential_temp()
    rk_stratiform_sedimentation()
    # RK2: detrain reserved liquid
    rk_stratiform_detrain_convective_condensate()
    apply_constituent_tendencies()
    qneg()
    geopotential_temp()
    # RK3a: cloud fractions (CAM4 version)
    convective_cloud_cover()
    convective_cloud_cover_diagnostics()
    compute_cloud_fraction()
    rk_stratiform_cloud_fraction_perturbation()
    rk_stratiform_cloud_fraction_perturbation_diagnostics()
    # RK3b: external forcings from non-macro/microphysics
    rk_stratiform_external_forcings()
    cloud_fraction_fice()               # ice fraction of stratus condensate/precipitation
    # RK3c1: prognostic cloud water, part one (pre-repartitioning values)
    prognostic_cloud_water()
    # RK3d1: stratiform condensate repartitioning
    rk_stratiform_condensate_repartioning()
    rk_stratiform_condensate_repartioning_diagnostics()
    apply_constituent_tendencies()
    qneg()
    geopotential_temp()
    # RK3c2 + RK3d2: prognostic cloud water, part two
    rk_stratiform_prognostic_cloud_water_tendencies()
    rk_stratiform_prognostic_cloud_water_tendencies_diagnostics()
    apply_constituent_tendencies()
    apply_heating_rate()
    qneg()
    geopotential_temp()


@ccpp_composite
def rk_stratiform():
    tropopause_find()
    rk_stratiform_diagnostics()             # initialise all RK diagnostics
    rk_stratiform_check_qtlcwat()           # QCWAT, TCWAT, LCWAT if not from snapshot
    energy_checked(process=_rk1_to_rk3,
                   flux_correction=rk_stratiform_prepare_flux_for_check_energy)
    # RK4: cloud fractions (third call to cldfrc)
    compute_cloud_fraction()
    compute_cloud_fraction_diagnostics()
    # RK5: cloud water and ice particle sizes for radiation
    rk_stratiform_cloud_optical_properties()
    rk_stratiform_cloud_optical_properties_diagnostics()
    # RK6: save non-micro and non-macrophysical external advective forcings
    rk_stratiform_save_qtlcwat()


rasch_kristjansson = Parameterization(name="Rasch-Kristjansson macro- and microphysics",
                                      run=rk_stratiform)
