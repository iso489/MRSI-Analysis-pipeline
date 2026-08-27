# MRSI Analysis Pipeline

Official Magnetic Resonance Spectroscopic Imaging (MRSI) analysis pipeline guide from the **Brain Tumor Research Group** at the **Montreal Neurological Institute, Montreal, Canada!**.

## Mission
We aim to **democratize access to multi-voxel MRSI** and support **its clinical integration**.

## What this is (and isn't)

The two halves of this repository have different maturity, and the distinction
matters if you intend to reproduce the paper.

**`Transcriptomic_Analysis/` is an executable pipeline.** Scripts 00-05 run end
to end from the two public Xena downloads; `verify_reproduction.R` checks a
fresh run against `EXPECTED_RESULTS.tsv` and exits non-zero on any mismatch.
Inputs are pinned by SHA-256 in `INPUT_MANIFEST.md`, package versions are
recorded in `dependencies.tsv`, and the complete sample-to-donor mapping is
published as `sample_donor_manifest.tsv`.

**`Spectroscopy/` is a documented procedure, not an executable pipeline.**
`PREPROCESSING.md` specifies the recommended steps and the FSL-MRS commands that
implement them; `mrsi_coregistration.py` is runnable and fully parameterised.
Reproducing the spectroscopic results additionally requires the vendor research
sequence and licensed processing software (MIDAS), which cannot be distributed
here. Do not expect to reproduce the 2-HG concentrations from this repository
alone.

## Inference is performed at the donor level

The 289 GTEx cortical samples used as the normal comparator are **not 289
independent brains**: they come from 164 donors, because GTEx sampled up to
three cortical regions per brain, and 91 donors contribute repeated
observations. Scripts 01-04 are sample-level and are retained for provenance.
**`05_donor_level_transcriptomics.R` performs the primary inference**, summing
each donor's counts and reporting two sensitivity analyses (one sample per donor
under a fixed region priority, and each cortical region separately).

## Reproducing the transcriptomic analysis

```bash
bash Transcriptomic_Analysis/download_xena.sh          # public inputs
shasum -a 256 -c Transcriptomic_Analysis/INPUT_SHA256  # verify what you got
sbatch Transcriptomic_Analysis/RNA_seq_analysis_script.sbatch
```

The Slurm script runs 00 -> 01 -> 02 -> 03 -> 04 -> 05 and then the verification
step. Set `GVN_BASE` to override the default `/scratch/$USER/glioma_vs_normal`.

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
| `00_prepare_pathways.R` | builds the 11 prespecified gene sets (five BRETIGEA neuronal-marker tiers, used unchanged from the published ranking, plus six KEGG/Reactome-based curated metabolic pathways) into `pathways.gmt`; the frozen gene lists used for the reported analyses are committed as `Transcriptomic_Analysis/pathways.gmt` |
| `01_deseq2_glioma_vs_normal.R` | gene filtering, protein-coding restriction (Ensembl release 100, from `ensembl100_gene_biotypes.tsv`) and DESeq2 `~ condition` differential expression; writes both `deseq2_results_glioma_vs_normal.tsv` and `deseq2_results_protein_coding.tsv` |
| `02_fgsea_analysis.R` | pre-ranked GSEA on the DESeq2 Wald statistic (`fgseaMultilevel`, `nPermSimple = 10000`) |
| `03_gsva_analysis.R` | per-sample pathway activity (GSVA, Gaussian kernel on variance-stabilized counts) |
| `04_negative_control_gsea.R` | calibrates each prespecified set against 500 size-matched random gene sets drawn from the same ranked universe |
| `05_donor_level_transcriptomics.R` | **primary inference**: within-donor aggregation (164 GTEx donors vs 662 TCGA participants), plus one-sample-per-donor and region-specific sensitivity analyses |
| `verify_reproduction.R` | checks a fresh run against `EXPECTED_RESULTS.tsv`; exits non-zero on mismatch |

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

