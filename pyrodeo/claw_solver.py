# -*- coding: utf-8 -*-
"""Defines a generic conservation law solver.
"""

from __future__ import print_function

import numpy as np

class ClawSolver:
    """Generic conservation law solver.

    Note:
        Serves as a base class to construct various solvers. Can not be used on its own.

    The following method is available:

    """
    def limiter(self, a, b, sb):
        """Limiter function to limit slopes/fluxes.

        This limiter function can, based on the parameter sb, emulate total variation diminishing limiters from minmod (sb = 1) to superbee (sb = 2).

        Args:
            a (ndarray): First slope to compare.
            b (ndarray): Second slope to compare.
            sb (float): Limiter parameter. Should be >= 1 (minmod, most diffusive limiter) and <= 2 (superbee, least diffusive limiter).

        Returns:
            ndarray: Limited slopes.

        """
        return (np.maximum(0.0*a,
                           np.minimum(sb*a,
                                      np.maximum(b, np.minimum(a, sb*b)))) +
                np.minimum(0.0*a,
                           np.maximum(sb*a,
                                      np.minimum(b, np.maximum(a, sb*b)))))
