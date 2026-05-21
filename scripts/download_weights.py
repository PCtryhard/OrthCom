#!/usr/bin/env python3
"""Download RFdiffusion model checkpoints into models/weights/."""

import sys
from pathlib import Path
from urllib.request import urlretrieve
from urllib.error import HTTPError


# Baker Lab download URLs for core RFdiffusion checkpoints
CHECKPOINTS = {
    "Base_ckpt.pt": "http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt",
    "Complex_base_ckpt.pt": "http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt",
}

DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent.parent / "models" / "weights"


def report_hook(block_num: int, block_size: int, total_size: int) -> None:
    """Print a simple progress bar for urlretrieve."""
    if total_size > 0:
        downloaded = block_num * block_size
        percent = min(downloaded * 100 / total_size, 100)
        print(f"\r  Progress: {percent:.1f}%", end="", flush=True)


def download_file(url: str, dest: Path) -> None:
    """Download a single file from *url* to *dest*."""
    print(f"Downloading {dest.name} ...")
    try:
        urlretrieve(url, dest, reporthook=report_hook)
        print()  # newline after progress bar
    except HTTPError as exc:
        print(f"\nERROR: Failed to download {dest.name} (HTTP {exc.code}).", file=sys.stderr)
        if dest.exists():
            dest.unlink()
        sys.exit(1)
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        if dest.exists():
            dest.unlink()
        sys.exit(1)


def main() -> None:
    output_dir = DEFAULT_OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    for filename, url in CHECKPOINTS.items():
        dest_path = output_dir / filename
        if dest_path.exists():
            print(f"Skipping {filename}: already exists at {dest_path}")
            continue
        download_file(url, dest_path)
        print(f"Saved to {dest_path}")

    print("All downloads complete.")


if __name__ == "__main__":
    main()
