# MRSI Analysis Pipeline (MNI-MRSI)

Official Magnetic Resonance Spectroscopic Imaging (MRSI) analysis pipeline guide from the **Brain Tumor Research Group** at the **Montreal Neurological Institute (MNI)**.

## Mission
We aim to **democratize access to multi-voxel MRSI** and support **its clinical integration**.

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

## Dependencies
- Python 3.8+
- FSL-MRS v2.4.9
- ANTs v2.6.2
- HD-BET
- FID-A (MATLAB)

## Recommended structure for your analysis

Tip: Keep your project reproducible by organizing your own analysis folder like this:
```text
project/
└── Subject_001/
    ├── Spectroscopy_Data/
    │   ├── Raw_Data/
    │   ├── Converted_Data/
    │   └── Processed_Data/
    │       ├── Frequency_and_Phase_Corrected_Data/
    │       ├── ECC_Data/
    │       └── Water_Removed_Data/
    ├── Anatomical_Data/
    │   ├── Raw_Data/
    │   ├── Converted_Data/
    │   ├── Denoised_Data/
    │   └── Brain_Extracted_Data/
    ├── PVE_Maps/
    └── Results/
```

## Citation
If you use this workflow in academic work, please cite:
- Oultache et al., 2025 [manuscript in preparation]

## License
This work is free to use (Apache-2.0 license)!

