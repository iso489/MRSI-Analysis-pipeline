# Input manifest

Every input this pipeline consumes, with the exact bytes it was run against.
Verify before running: `shasum -a 256 -c INPUT_SHA256`.

## Downloaded inputs (UCSC Xena / Toil recompute hub)

| File | Bytes | SHA-256 |
|---|---|---|
| `TcgaTargetGtex_gene_expected_count.gz` | 1,289,124,490 | `ee706cb15448417dc522f1efb2dbdd34ef342c0d5a038cb93ecbefa7990903db` |
| `TcgaTargetGTEX_phenotype.txt.gz` | 135,753 | `ba4d4461cff0fe5e91cc3e58f793aa47ff3c5d6fbf77d65ccff5a07231aaf2db` |

Base URL: `https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/`
Both files carry `Last-Modified: 2021-04-09`; the hashes above are what the
published analysis was run against.

**Value encoding.** The expression matrix stores `log2(RSEM expected count + 1)`
rounded to four decimal places. Script 01 back-transforms with
`round(2^x - 1)`. That reconstruction is exact for 97.7% of matrix entries; the
remaining 2.3% admit more than one integer, and the worst-case relative error
over the whole matrix is 3.5e-5. DESeq2 dispersions are orders of magnitude
larger, so the rounding is immaterial, but it is stated rather than assumed.

## Derived input you must build

`ensembl100_gene_biotypes.tsv` — two columns, `gene_id` and `biotype`, from the
Ensembl release-100 GTF. Script 01 requires it and **fails loudly** if absent.

```bash
curl -sL https://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz \
| gunzip -c | awk -F'\t' '$3=="gene"' \
| sed -n 's/.*gene_id "\([^"]*\)".*gene_biotype "\([^"]*\)".*/\1\t\2/p' \
| sort -u | sed '1i gene_id\tbiotype' > ensembl100_gene_biotypes.tsv
```

Note that Xena gene IDs carry a version suffix (`ENSG00000000003.14`). Match on
the unversioned identifier if your GTF-derived table lacks the suffix.

## Cohort selection

| Group | Selection | n |
|---|---|---|
| Glioma | TCGA GBM + LGG, `_sample_type == "Primary Tumor"` | 662 samples from 662 participants |
| Normal | GTEx `Brain - Frontal Cortex (Ba9)`, `Brain - Anterior Cingulate Cortex (Ba24)`, `Brain - Cortex` | 289 samples from **164 donors** |

The 290th eligible phenotype record, `GTEX-13S7M-0011-R10b-SM-5PNZB` (BA9), has
no column in the expression matrix, which is why 289 and not 290 samples are
analysed.

**The 289 GTEx samples are not 289 independent brains.** GTEx sampled up to
three cortical regions per donor: 73 donors contributed one sample, 59 two, 31
three and one (`GTEX-NPJ8`) five, so 91 donors contribute repeated observations.
`sample_donor_manifest.tsv` gives the complete sample -> donor -> region mapping.
Script `05_donor_level_transcriptomics.R` performs the primary inference at the
donor level; scripts 01-04 remain sample-level and are retained for provenance.
