"""Category: drag from the surface that is not part of vertical diffusion."""

from xdsl_ccpp.frontend.py_api import ccpp_composite

from categories.interface import run_setup


@ccpp_composite
def surface_drag(*, form_drag):
    """Roles:
        form_drag: orographic form drag from subgrid orography (required).
    """
    run_setup(form_drag)
    form_drag.run()
