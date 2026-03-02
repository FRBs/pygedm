#!/usr/bin/env python3
"""Download gedm_dist_maps.hkl from Zenodo if not already present."""

import sys
import os
from pathlib import Path
import urllib.request
import urllib.error

# Configuration
DATA_DIR = Path("/app/assets")
DATA_FILE = DATA_DIR / "gedm_dist_maps.hkl"
ZENODO_URL = "https://zenodo.org/records/18779007/files/gedm_dist_maps.hkl?download=1"
MIN_SIZE = 1024 * 1024  # 1 MB — file should be ~60 MB, so anything smaller is corrupt


def main():
    # Check for a local copy in the build context (mounted at /tmp/ctx)
    local_path = Path("/tmp/ctx/data/gedm_dist_maps.hkl")
    if local_path.exists() and local_path.stat().st_size > MIN_SIZE:
        print(f"✓ Found local copy in build context ({local_path.stat().st_size / 1024**2:.1f} MB)")
        print(f"  Copying to {DATA_FILE}...")
        import shutil
        shutil.copy2(local_path, DATA_FILE)
        return 0

    # Check if file already exists in container
    if DATA_FILE.exists() and DATA_FILE.stat().st_size > MIN_SIZE:
        print(f"✓ {DATA_FILE} already exists ({DATA_FILE.stat().st_size / 1024**2:.1f} MB)")
        return 0

    # File is missing or too small — download from Zenodo
    print(f"↓ Downloading gedm_dist_maps.hkl from Zenodo (~60 MB)...")
    print(f"  URL: {ZENODO_URL}")

    try:
        urllib.request.urlretrieve(ZENODO_URL, DATA_FILE)
    except urllib.error.URLError as e:
        print(f"✗ FATAL: Failed to download from Zenodo: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"✗ FATAL: Unexpected error during download: {e}", file=sys.stderr)
        return 1

    # Verify the downloaded file
    if not DATA_FILE.exists():
        print(f"✗ FATAL: Download appeared to succeed but file does not exist", file=sys.stderr)
        return 1

    size = DATA_FILE.stat().st_size
    if size < MIN_SIZE:
        print(f"✗ FATAL: Downloaded file is only {size / 1024**2:.1f} MB (expected ~60 MB)", file=sys.stderr)
        return 1

    # Verify the file is readable HDF5
    try:
        import h5py
        with h5py.File(DATA_FILE, 'r') as h:
            print(f"✓ Downloaded successfully ({size / 1024**2:.1f} MB) and verified as valid HDF5")
    except ImportError:
        print(f"✓ Downloaded successfully ({size / 1024**2:.1f} MB)")
    except Exception as e:
        print(f"✗ FATAL: File is not valid HDF5: {e}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
