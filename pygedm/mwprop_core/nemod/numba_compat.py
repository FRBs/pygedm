# mwprop v2.0 Jan 2026

"""
Numba compatibility module

Provides a fallback njit decorator when numba is unavailable.
This allows code to work without numba installed, though without
the performance benefits of JIT compilation.
"""

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    def njit(*args, **kwargs):
        """Dummy decorator when numba unavailable"""
        def decorator(func):
            return func
        # Handle @njit (called with function) vs @njit() (called with args)
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]
        return decorator

__all__ = ['njit', 'HAS_NUMBA']
