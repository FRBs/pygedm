# Migration from pkg_resources to importlib

## Summary

This package has been updated to use modern Python packaging standards, replacing the deprecated `pkg_resources` API with `importlib.metadata` and `importlib.resources`.

## Changes Made

### 1. Version Detection
**Before:**
```python
from pkg_resources import get_distribution, DistributionNotFound
__version__ = get_distribution("pygedm").version
```

**After:**
```python
from importlib.metadata import version, PackageNotFoundError
__version__ = version("pygedm")
```

### 2. Resource File Access
**Before:**
```python
from pkg_resources import resource_filename
DATAPATH = os.path.dirname(resource_filename("pygedm", "spiral.txt"))
```

**After:**
```python
from importlib.resources import files
DATAPATH = str(files("pygedm"))
```

### 3. Updated Dependencies
- Added `importlib-metadata` for Python < 3.8 (backport)
- Added `importlib-resources` for Python < 3.9 (backport)
- Minimum Python version remains 3.8+

### 4. Modern Packaging
- Added `pyproject.toml` with PEP 621 metadata
- Updated `setup.py` with new dependencies
- Updated `requirements.txt`

## Benefits

1. **No more pkg_resources errors**: Works with setuptools 81.0+ and modern Python environments
2. **Standard library approach**: Uses built-in modules for Python 3.8+
3. **Future-proof**: Follows current Python packaging best practices
4. **Faster imports**: `importlib.metadata` is more efficient than `pkg_resources`

## Installation

The package now requires Python 3.8 or later. Install as usual:

```bash
pip install pygedm
```

For development:
```bash
pip install -e .
```

## Notes

- The `pkg_resources` API was deprecated and is being removed from setuptools
- All functionality remains the same; only internal implementation changed
- No API changes for end users
