"""Held-Suarez (1994) suite: Python version of suites/suite_held_suarez_1994.xml.

Dry dynamical-core test with idealized forcing.
"""

from xdsl_ccpp.frontend.py_api import ccpp_suite, export_suite

from categories.diagnostics import state_diagnostics, tendency_diagnostics
from categories.idealized import idealized_physics
from parameterizations.held_suarez import held_suarez


# In a CCPP suite class, each method is a group recorded by @ccpp_suite,
# not a real method: it takes no 'self', and the class holds only groups.
# pylint: disable=no-method-argument,too-few-public-methods
@ccpp_suite("held_suarez_1994", version="1.0")
class HeldSuarez1994Suite:
    """The 'held_suarez_1994' suite."""

    def physics_before_coupler():
        """Physics before the coupler."""
        idealized_physics(scheme=held_suarez)
        state_diagnostics()

    def physics_after_coupler():
        """Physics after the coupler."""
        tendency_diagnostics()


if __name__ == "__main__":
    export_suite(HeldSuarez1994Suite)
