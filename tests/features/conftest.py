"""Shared world for the acceptance scenarios.

The scenarios are written for a curator, so the steps carry the whole burden
of translating clinical language into the system's vocabulary. Everything the
steps need to remember between one step and the next lives in ``mundo``.
"""

from __future__ import annotations

import pytest


class Mundo(dict):
    """State passed between steps of one scenario."""


@pytest.fixture()
def mundo() -> Mundo:
    return Mundo()
