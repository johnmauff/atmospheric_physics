"""RRTMGP test suite: Python version of test/test_suites/suite_rrtmgp.xml.

Radiation only, for testing RRTMGP on its own: the test scheme
initialize_constituents sets the constituents a standalone radiation test
needs, then RRTMGP runs with zero aerosol optics.  Used by CAM-SIMA's
composite RRTMGP tests (CCPP_SUITE_DIR=suites-py).
"""

from xdsl_ccpp.frontend.py_api import ccpp_scheme_ref, ccpp_suite, export_suite

from categories.radiation import radiation
from parameterizations.aerosol_optics import zero_aerosol_optics
from parameterizations.rrtmgp import rrtmgp_shortwave_aerosols

initialize_constituents = ccpp_scheme_ref("initialize_constituents")   # test scheme


# In a CCPP suite class, each method is a group recorded by @ccpp_suite,
# not a real method: it takes no 'self', and the class holds only groups.
# pylint: disable=no-method-argument,too-few-public-methods
@ccpp_suite("rrtmgp", version="1.0")
class RrtmgpSuite:
    """The 'rrtmgp' suite."""

    def physics_after_coupler():
        """Physics after the coupler: radiation only."""
        initialize_constituents()
        radiation(scheme=rrtmgp_shortwave_aerosols, aerosol_optics=zero_aerosol_optics)


if __name__ == "__main__":
    export_suite(RrtmgpSuite)
