"""Category: atmospheric chemistry."""

from xdsl_ccpp.frontend.py_api import ccpp_composite

from categories.interface import run_setup


@ccpp_composite
def chemistry(*, scheme):
    """Roles:
        scheme: the chemistry package (required).
    """
    run_setup(scheme)
    scheme.run()
