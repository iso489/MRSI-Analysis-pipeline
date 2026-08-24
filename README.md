# MRSI Analysis Pipeline

Official Magnetic Resonance Spectroscopic Imaging (MRSI) analysis pipeline guide from the **Brain Tumor Research Group** at the **Montreal Neurological Institute, Montreal, Canada!**.

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

## Transcriptomic analysis

`Transcriptomic_Analysis/` holds the R scripts behind the radiogenomic context
analysis (TCGA glioma versus GTEx normal cortex, UCSC Xena Toil recompute hub):

| script | what it does |
|---|---|
| `00-Prepare_Pathways.R` | builds the 11 prespecified gene sets (five BRETIGEA neuronal-marker tiers, used unchanged from the published ranking, plus six KEGG/Reactome-based curated metabolic pathways) into `pathways.gmt`; the frozen gene lists used for the reported analyses are committed as `Transcriptomic_Analysis/pathways.gmt` |
| `01_deseq2_analysis.R` | gene filtering, protein-coding restriction (Ensembl release 100) and DESeq2 `~ condition` differential expression |
| `02_fgsea_analysis.R` | pre-ranked GSEA on the DESeq2 Wald statistic (`fgseaMultilevel`, `nPermSimple = 10000`) |
| `03_gsva_analysis.R` | per-sample pathway activity (GSVA, Gaussian kernel on variance-stabilized counts) |
| `04_negative_control_gsea.R` | calibrates each prespecified set against 500 size-matched random gene sets drawn from the same ranked universe |

Script 04 exists because 92% of tested protein-coding genes are differentially
expressed between bulk glioma and normal cortex. With a genome-wide shift that
large, any gene set can look enriched, so the calibration tests that objection
empirically rather than arguing about it. It is the analysis the manuscript
relies on when it states which pathway results are interpretable.

`RNA_seq_analysis_script.sbatch` is the Slurm job used to run the pipeline on
the Digital Research Alliance of Canada infrastructure.

## Citation
If you use this workflow in academic work, please cite the relevant paper:
- [manuscript in preparation]

## License
This work is free to use (Apache-2.0 license)!

