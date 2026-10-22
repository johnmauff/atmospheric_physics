"""MUSICA suite: Python version of suites/suite_musica.xml.

Chemistry only.
"""

from xdsl_ccpp.frontend.py_api import ccpp_suite, export_suite

from categories.chemistry import chemistry
from categories.diagnostics import state_diagnostics
from parameterizations.musica import musica


# In a CCPP suite class, each method is a group recorded by @ccpp_suite,
# not a real method: it takes no 'self', and the class holds only groups.
# pylint: disable=no-method-argument,too-few-public-methods
@ccpp_suite("musica", version="1.0")
class MusicaSuite:
    """The 'musica' suite."""

    def physics_after_coupler():
        """Physics after the coupler."""
        chemistry(scheme=musica)
        state_diagnostics()


if __name__ == "__main__":
    export_suite(MusicaSuite)
