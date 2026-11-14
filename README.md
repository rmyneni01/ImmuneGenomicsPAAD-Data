# FGA and immune checkpoint / cytolytic gene expression in PAAD (TCGA)

Code for: Myneni, R. *Characterizing the Association of Genomic Instability with Immune Checkpoint and Cytolytic Gene Expression in Pancreatic Adenocarcinoma* (2026).

## Data
Public, de-identified TCGA data from cBioPortal, study `paad_tcga_pan_can_atlas_2018` (n = 184; RNA-seq for 177).
- `data/paad_tcga_pan_can_atlas_2018_clinical_data.tsv`: clinical + FGA, mutation count, TMB
- `data/mRNA_Expression__RSEM__Batch_normalized_...txt`: RSEM for 8 genes (CD274, PDCD1, CTLA4, CD8A, GZMB, PRF1, LAG3, HAVCR2)
- ABSOLUTE purity: `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` (GDC PanCanAtlas publications page, Taylor et al. 2018)
- Methylation-based leukocyte fraction: `TCGA_all_leuk_estimate.masked.20170107.tsv` (GDC PanImmune page, Thorsson et al. 2018)

## Run
```bash
pip install pandas numpy scipy lifelines
python analysis.py \
  --clinical data/paad_tcga_pan_can_atlas_2018_clinical_data.tsv \
  --expr data/mRNA_Expression__RSEM__Batch_normalized_from_Illumina_HiSeq_RNASeqV2___1_.txt \
  --purity data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt --purity-id-col sample --purity-col purity \
  --leuk data/TCGA_all_leuk_estimate.masked.20170107.tsv   # both optional
```
Outputs `results/correlations.csv`, `results/fga_vs_leukocyte_fraction.csv` (Spearman and partial Spearman, BH-FDR within each analysis) and `results/cox_cyt.csv`.

## Analyses
- FGA vs expression: all tumors, PDAC-annotated, FGA thresholds (0 to 0.10)
- Partial Spearman controlling for mutation count and, if supplied, purity
- Cox model: log2 CYT + age + sex + stage (+ purity)
- CYT = geometric mean of GZMB and PRF1

## Notes
- Expression file from cBioPortal may have lost line breaks; the script handles this.
- Without `--purity`, purity-adjusted analyses are skipped.
- Data: TCGA Research Network (https://www.cancer.gov/tcga).
