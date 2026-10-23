"""Kessler test suite: Python version of test/test_suites/suite_kessler_test.xml.

Kessler microphysics on its own, with state diagnostics: no energy check
and no physics after the coupler. Used by CAM-SIMA's composite Kessler tests
(CCPP_SUITE_DIR=test/test_suites-py).
"""

from xdsl_ccpp.frontend.py_api import ccpp_suite, export_suite

from categories.diagnostics import state_diagnostics
from parameterizations.kessler import kessler


# In a CCPP suite class, each method is a group recorded by @ccpp_suite,
# not a real method: it takes no 'self', and the class holds only groups.
# pylint: disable=no-method-argument,too-few-public-methods
@ccpp_suite("kessler_test", version="1.0")
class KesslerTestSuite:
    """The 'kessler_test' suite."""

    def physics_before_coupler():
        """Physics before the coupler: Kessler, then state diagnostics."""
        kessler.run()
        state_diagnostics(parameterizations=[kessler])


if __name__ == "__main__":
    export_suite(KesslerTestSuite)
