# FGA and immune checkpoint / cytolytic gene expression in PAAD (TCGA)

Code for: Myneni, R. *Characterizing the Association of Genomic Instability with Immune Checkpoint and Cytolytic Gene Expression in Pancreatic Adenocarcinoma* (2026).

## Data
Public, de-identified TCGA data from cBioPortal, study `paad_tcga_pan_can_atlas_2018` (n = 184; RNA-seq for 177).
- `data/paad_tcga_pan_can_atlas_2018_clinical_data.tsv`: clinical + FGA, mutation count, TMB
- `data/mRNA_Expression__RSEM__Batch_normalized_...txt`: RSEM for 8 genes (CD274, PDCD1, CTLA4, CD8A, GZMB, PRF1, LAG3, HAVCR2)
- Tumor purity (not in cBioPortal clinical file): Aran et al. 2015, *Nat Commun* 6:8971, Supplementary Data (CPE column), keyed by TCGA sample barcode.

## Run
```bash
pip install pandas numpy scipy lifelines
python analysis.py \
  --clinical data/paad_tcga_pan_can_atlas_2018_clinical_data.tsv \
  --expr data/mRNA_Expression__RSEM__Batch_normalized_from_Illumina_HiSeq_RNASeqV2___1_.txt \
  --purity data/purity.csv --purity-col CPE   # optional
```
Outputs `results/correlations.csv` (Spearman and partial Spearman, BH-FDR within each analysis) and `results/cox_cyt.csv`.

## Analyses
- FGA vs expression: all tumors, PDAC-annotated, FGA thresholds (0 to 0.10)
- Partial Spearman controlling for mutation count and, if supplied, purity
- Cox model: log2 CYT + age + sex + stage (+ purity)
- CYT = geometric mean of GZMB and PRF1

## Notes
- Expression file from cBioPortal may have lost line breaks; the script handles this.
- Without `--purity`, purity-adjusted analyses are skipped.
- Data: TCGA Research Network (https://www.cancer.gov/tcga).
