# mwprop.nemod v2.0 Feb 2026

"""
galaxy_model.py
===============
GalaxyModel — container for pre-computed, call-invariant NE20x model data.

Usage
-----
The default singleton is used automatically by ``ne_arms_ne2001p``:

    from pygedm.mwprop_core.nemod.density import density_2001_smooth_comps  # works as before

To switch NE2001 ↔ NE2025:

    from pygedm.mwprop_core.nemod.config_nemod import set_model
    set_model('ne2025')        # updates globals AND refreshes default_model

To create a custom instance (e.g. for testing a modified parameter set):

    from pygedm.mwprop_core.nemod.galaxy_model import GalaxyModel
    my_model = GalaxyModel()

A module-level singleton ``default_model`` is created on first import.
"""

import numpy as np
from scipy.interpolate import CubicSpline


class GalaxyModel:
    """
    Pre-computes all call-invariant data for the NE20x galaxy model.

    The heavy construction is done once in ``__init__``; every subsequent
    call to ``ne_arms_ne2001p`` performs only cheap attribute lookups.

    Attributes
    ----------
    narms : int
        Number of spiral arms.
    Ncoarse : int
        Number of coarse angular samples used to define each arm.
    arm_radius_splines : list[CubicSpline]
        ``CubicSpline(theta, radius)`` for each arm — built once, reused
        across all LoS evaluations instead of being recreated per call.
    arm_index : np.ndarray, shape (narms,), int
        Maps sequential arm index *j* to the TC93 / NE2001 arm numbering.
    warm, harm, narm : np.ndarray, shape (narms,)
        Per-arm width, height and amplitude scale factors from *Dgal*.
    wa, Aa, ha, Fa, na : float
        Scalar arm parameters from *Dgal*.
    """

    def __init__(self):
        self._setup()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _setup(self):
        """Extract / build all pre-computable model data from config globals."""
        # Late import avoids a circular dependency at module-load time:
        #   galaxy_model → config_nemod  is fine
        #   config_nemod → galaxy_model  must NOT happen at module level
        from pygedm.mwprop_core.nemod.config_nemod import (
            Dgal, Darmmap, coarse_arms, armsplines, th1, r1,
        )

        narms   = coarse_arms.shape[1]
        Ncoarse = coarse_arms.shape[2]
        self.narms   = narms
        self.Ncoarse = Ncoarse

        # ----------------------------------------------------------------
        # Spiral arm radius splines
        # ``config_nemod.setup_spiral_arms()`` already builds ``armsplines``.
        # We simply expose them as a plain list so ``ne_arms_ne2001p`` can
        # drop the ``if 'armsplines' in globals()`` guard entirely.
        # ----------------------------------------------------------------
        if armsplines is not None and len(armsplines) == narms:
            self.arm_radius_splines = list(armsplines)
        else:
            # Fallback in case setup_spiral_arms was not yet called.
            self.arm_radius_splines = [
                CubicSpline(th1[j, :], r1[j, :]) for j in range(narms)
            ]

        # ----------------------------------------------------------------
        # Integer arm-index mapping
        # Previously: np.fromiter(...) executed on every ne_arms call.
        # ----------------------------------------------------------------
        self.arm_index = np.fromiter(
            (Darmmap[str(j)] for j in range(narms)), dtype=int
        )

        # ----------------------------------------------------------------
        # Per-arm scale-parameter arrays
        # Previously: three list-comprehension+dict-lookup passes executed
        # inside the per-arm loop on every ne_arms call.
        # ----------------------------------------------------------------
        self.warm = np.array([Dgal['warm' + str(jj)] for jj in self.arm_index])
        self.harm = np.array([Dgal['harm' + str(jj)] for jj in self.arm_index])
        self.narm = np.array([Dgal['narm' + str(jj)] for jj in self.arm_index])

        # ----------------------------------------------------------------
        # Scalar arm parameters
        # Previously: four dict lookups executed on every ne_arms call.
        # ----------------------------------------------------------------
        self.wa = Dgal['wa']
        self.Aa = Dgal['Aa']
        self.ha = Dgal['ha']
        self.Fa = Dgal['Fa']
        self.na = Dgal['na']

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def refresh(self):
        """Re-run setup; call after ``set_model()`` switches the galaxy model.

        ``config_nemod.set_model()`` calls this automatically via a lazy
        import so callers normally never need to invoke it directly.
        """
        self._setup()

    def __repr__(self):
        return (
            f"GalaxyModel(narms={self.narms}, Ncoarse={self.Ncoarse}, "
            f"wa={self.wa:.4f}, Aa={self.Aa:.4f})"
        )


# ---------------------------------------------------------------------------
# Module-level singleton — created once when this module is first imported.
# ``ne_arms_ne2001p`` uses this by default; ``set_model()`` refreshes it.
# ---------------------------------------------------------------------------
default_model = GalaxyModel()
