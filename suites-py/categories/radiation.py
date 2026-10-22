"""Category: radiation."""

from xdsl_ccpp.frontend.py_api import ccpp_composite

from categories.interface import nothing, run_setup


@ccpp_composite
def radiation(*, scheme, cloud_water=None, aerosol_optics=None):
    """Radiative transfer and the heating it produces.

    Roles:
        scheme:         the radiation scheme (required), a Parameterization
                        whose run takes ``aerosol_optics=`` -- the scheme
                        decides where in its sequence aerosol optics run.
        cloud_water:    cloud water paths radiation reads, computed first;
                        or None when the cloud scheme provides them.
        aerosol_optics: the model's aerosol optical properties, or None.
    """
    if cloud_water is not None:
        run_setup(cloud_water)
        cloud_water.run()
    run_setup(scheme)
    scheme.run(aerosol_optics=aerosol_optics.run if aerosol_optics is not None else nothing)
