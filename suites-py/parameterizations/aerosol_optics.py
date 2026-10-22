"""Aerosol optical properties for radiation."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization

aerosol_optics_scheme      = ccpp_scheme_ref("aerosol_optics")
aerosol_optics_diagnostics = ccpp_scheme_ref("aerosol_optics_diagnostics")
compute_relative_humidity  = ccpp_scheme_ref("compute_relative_humidity")
aerosol_optics_zero_stub   = ccpp_scheme_ref("aerosol_optics_zero_stub")


@ccpp_composite
def bulk_aerosol_optical_properties():
    """Bulk aerosols.  The relative humidity they need comes from the cloud
    fraction parameterization, which must run first."""
    aerosol_optics_scheme()
    aerosol_optics_diagnostics()


@ccpp_composite
def bulk_aerosol_optical_properties_with_humidity():
    """Bulk aerosols, computing the relative humidity they need first
    (CAM7: no cloud scheme provides it)."""
    compute_relative_humidity()
    aerosol_optics_scheme()
    aerosol_optics_diagnostics()


bulk_aerosol_optics = Parameterization(name="bulk aerosol optics",
                                       run=bulk_aerosol_optical_properties)
# No aerosols: zero aerosol optical properties (RRTMGP test suite).
zero_aerosol_optics = Parameterization(name="no aerosols (zero optics)",
                                       run=aerosol_optics_zero_stub)
bulk_aerosol_optics_with_humidity = Parameterization(
    name="bulk aerosol optics, with relative humidity",
    run=bulk_aerosol_optical_properties_with_humidity)
