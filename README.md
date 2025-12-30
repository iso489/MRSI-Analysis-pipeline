# MRSI Analysis Pipeline (MNI)

Official Magnetic Resonance Spectroscopic Imaging (MRSI) analysis pipeline guide from the **Brain Tumor Research Group** at the **Montreal Neurological Institute (MNI)**.

## Mission
We aim to **democratize access to multi-voxel MRSI** and support **more standardized clinical integration**.

## What this is (and isn’t)
**This repository provides**: a documented, step-by-step pipeline (recommended tools, parameters, and workflow order) for multi-voxel MRSI analysis.
**This repository does not provide**: a one-click automated workflow (no single script that runs everything end-to-end).

## Pipeline overview
This guide describes a reproducible workflow typically including:
1. MRI preprocessing (e.g., denoising, brain extraction)
2. MRSI preprocessing (format conversion, frequency/phase correction, water removal, etc.)
3. Coregistration of MRSI to anatomical MRI
4. Basis set generation (sequence-matched)
5. Tissue segmentation / partial-volume estimation (GM/WM/CSF)
6. Spectral fitting and quality control (QC)
7. Reporting: metabolite maps and ROI summaries (depending on your use case)

## Tools referenced
Our workflow is based on open-source tools:
- **HD-BET** (brain extraction)
- **ANTs** (MRI preprocessing)
- **FSL-MRS** (MRS preprocessing and fitting)
- **FID-A** (basis set generation)

## How to use this repository
- Start here: `docs/` (or `pipeline.md`)  
- Follow the steps in order and adapt commands/parameters to your data format and scanner protocol.
- Use the included examples as templates (configs, command snippets, checklists).

## Recommended structure for your analysis
> Tip: Keep your project reproducible by organizing your own analysis folder like this:
  >project/
  >raw/
  >derivatives/
  >anat/
  >mrsi/
  >qc/
  >results/
  >notes/

## Citation
If you use this workflow in academic work, please cite:
- TODO: paper / DOI / preprint

## License
