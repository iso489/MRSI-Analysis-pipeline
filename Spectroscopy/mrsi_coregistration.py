#!/usr/bin/env python3
"""
Warp brain mask into MRSI (CSI) space.
Uses fslpy approach to create a dummy image matching CSI space exactly.
Relies on header-based alignment (applywarp --usesqform).

"""

import argparse
import sys
import numpy as np
from pathlib import Path
from fsl.data.image import Image
from fsl.wrappers import applywarp

# -----------------------------------------------------------------------------
# 1) Input/output paths, taken from the command line.
#
# These were previously literal "enter the path" placeholders, which made the
# script unrunnable as published. Every path is now an argument, so the script
# can be executed directly and driven from a batch loop.
# -----------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Warp a brain mask into MRSI (CSI) space using header-based "
                    "alignment (applywarp --usesqform).")
    ap.add_argument("--mrsi", type=Path, required=True,
                    help="CSI-space reference image (defines the target grid).")
    ap.add_argument("--brain-mask", type=Path, required=True,
                    help="Brain mask in structural space, to be warped into CSI space.")
    ap.add_argument("--output", type=Path, required=True,
                    help="Output directory; created if absent.")
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="Binarisation threshold for the warped mask (default: 0.5).")
    a = ap.parse_args(argv)
    for label, path in (("--mrsi", a.mrsi), ("--brain-mask", a.brain_mask)):
        if not path.exists():
            ap.error(f"{label} does not exist: {path}")
    a.output.mkdir(parents=True, exist_ok=True)
    a.output_mask = a.output / "brain_mask_in_csi_space.nii.gz"
    a.output_mask_bin = a.output / "brain_mask_in_csi_space_bin.nii.gz"
    return a

args = parse_args()

# -----------------------------------------------------------------------------
# 2) Make dummy nifti - create zeros matching CSI space exactly
# -----------------------------------------------------------------------------
print(f"Loading CSI reference: {args.mrsi}")
mrsi_in = Image(str(args.mrsi))

# Create dummy image as nothing (zeros) matching CSI dimensions
# Note: only use spatial dimensions [0:3], ignore any 4th dimension (time/spectral)
tmp_img = np.zeros(mrsi_in.shape[0:3])
tmp_img = Image(tmp_img, xform=mrsi_in.voxToWorldMat)

print(f"  CSI shape: {mrsi_in.shape[0:3]}")
print(f"  Dummy image created with matching geometry")

# -----------------------------------------------------------------------------
# 3) Warp brain mask into CSI space using the dummy as reference
# -----------------------------------------------------------------------------
print(f"\nWarping brain mask into CSI space...")
print(f"  Input mask: {args.brain_mask}")

if not args.brain_mask.exists():
    print(f"ERROR: Brain mask not found at {args.brain_mask}")
    exit(1)

# Use nearest neighbor interpolation to preserve binary mask values
applywarp(
    str(args.brain_mask),   # input mask to warp
    tmp_img,                # reference dummy image (CSI grid)
    str(args.output_mask),  # output
    usesqform=True,         # use sform/qform for alignment (no warp field)
    super=True,             # supersampling for better interpolation
    superlevel='a',
    interp='nn'             # nearest neighbor for binary mask
)

print(f"Brain mask warped to: {args.output_mask}")

# -----------------------------------------------------------------------------
# 4) Binarize the output mask to ensure clean 0/1 values
# -----------------------------------------------------------------------------
print(f"\nBinarizing the warped mask...")
mask_img = Image(str(args.output_mask))
mask_data = mask_img.data

# Binarise at the requested threshold (default 0.5)
mask_binary = (mask_data > args.threshold).astype(np.float32)
mask_binary_img = Image(mask_binary, xform=mask_img.voxToWorldMat)
mask_binary_img.save(str(args.output_mask_bin))

print(f"Binary mask saved to: {args.output_mask_bin}")

# -----------------------------------------------------------------------------
# 5) Summary
# -----------------------------------------------------------------------------
print("\n" + "="*70)
print("REGISTRATION COMPLETE")
print("="*70)
print(f"Output files:")
print(f"Original warped mask: {args.output_mask}")
print(f"Binarized mask:       {args.output_mask_bin}")
