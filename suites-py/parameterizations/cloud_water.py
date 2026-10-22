"""Cloud water paths for radiation (CAM5, CAM7): convective and total."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization

convective_cloud_water             = ccpp_scheme_ref("convective_cloud_water")
convective_cloud_water_diagnostics = ccpp_scheme_ref("convective_cloud_water_diagnostics")
cloud_water_paths_scheme           = ccpp_scheme_ref("cloud_water_paths")
cloud_water_paths_diagnostics      = ccpp_scheme_ref("cloud_water_paths_diagnostics")


@ccpp_composite
def cloud_water_for_radiation():
    convective_cloud_water()
    convective_cloud_water_diagnostics()
    cloud_water_paths_scheme()
    cloud_water_paths_diagnostics()


cloud_water_paths = Parameterization(name="cloud water paths", run=cloud_water_for_radiation)
