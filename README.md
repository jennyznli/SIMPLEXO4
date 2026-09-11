# SIMPLEXO4

## 1. select_variants_annotate.smk

Starting with Penn Medicine Biobank's (PMBB) pre-filtered GoldiLocks (GL) per-chromosome VCF, apply a few additional QC filters on top of existing filters (see PMBB Release 4.0 Whole-Exome Sequencing documentation):

- split/normalized multiallelics

- genotype-level no-calling (DP\<7/GQ\<20 -\> `./.`)

- SNP/indel-specific Allele Balance site removal

- exclusion of gender-discordant sample

Performs variant annotation and compiles a genome-wide CSV file with annotations.

#### Config/Inputs

+---------------------+---------------------------------------------------------+---------------+
|                     |                                                         |               |
+=====================+=========================================================+===============+
| **Key**             | **Description**                                         | **Required?** |
+---------------------+---------------------------------------------------------+---------------+
| `input.chromosomes` | Per-chromosome file path template (`{CHR}` placeholder) | yes           |
+---------------------+---------------------------------------------------------+---------------+
| `input.samples`     | Input file path w/ all case + control IDs, one per line | yes           |
+---------------------+---------------------------------------------------------+---------------+
| `input.sex_info`    | No header FID/IID/SEX file (1=male, 2=female)           | yes           |
+---------------------+---------------------------------------------------------+---------------+
| `input.targets`     | Optional BED/region file to restrict sites              | no            |
+---------------------+---------------------------------------------------------+---------------+

#### **Rules**

- `bcftools_filter_per_chr`
  - [input]{.underline}: per-chr GL VCF (config)
  - [output]{.underline}: data/bcftools/chr{CHR}.site-qc.bcf
  - `bcftools view -S -a` restricts to the sample set; `-a` trims ALT alleles/INFO not seen in remaining samples; optionally restricts to target sites via `-R`
  - `bcftools annotate --set-id` overwrites variant ID with CHROM_POS_REF_ALT
  - `bcftools +fill-tags` recomputes tags on this sample subset -- differs from PMBB's full-cohort values
  - `bcftools filter -s LowCallRate` filters out `F_MISSING > 0.05`; recomputed on this subset, not PMBB's full cohort
  - `bcftools view -f PASS` restricts to PASS sites
- `bcftools_no_sample_vcf`
  - [input]{.underline}: data/bcftools/chr{CHR}.site-qc.bcf
  - [output]{.underline}: data/bcftools/chr{CHR}.site-qc.no_sample.vcf
  - removes genotypes for a faster VEP pass
- `first_variant_annotation`
  - [input]{.underline}: data/bcftools/chr{CHR}.site-qc.no_sample.vcf
  - [output]{.underline}: data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf
  - runs VEP annotation in a singularity container and all plugins
- `plink2_annotation_pgen`
  - [input]{.underline}: data/bcftools/chr{CHR}.site-qc.bcf
  - [output]{.underline}: data/preprocess/chr{CHR}.annotation.{pgen,pvar,psam}
  - converts site-qc BCF to PLINK w/ double ID and updated sex
- `parse_annotation_no_sample_vep`
  - [input]{.underline}: data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf
  - [output]{.underline}: data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv
  - parses the VEP-annotated VCF into a CSV in `no_sample` mode
- `combine_vep_reports`
  - [input]{.underline}: data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv
  - [output]{.underline}: data/preprocess/all_chr.vep.report.csv
  - combines all per-chromosome VEP report CSVs into one file (input to `preprocess_regenie.py`)

#### **Submission**

``` bash
bsub -N -J select_variants_annotate -eo logs/select_variants_annotate.e -oo logs/select_variants_annotate.o \
    snakemake -s scripts/select_variants_annotate.smk --configfile config/select_variants.yaml --workflow-profile \
    profiles/lsf -j 22 
```

#### **Outputs**

- data/preprocess/chr{CHR}.annotation.{pgen,pvar,psam}
- data/preprocess/all_chr.vep.report.csv

## 2. regenie.smk

Full Regenie pipeline: builds the sample list, prepares Step 1 genotypes (array or LD-pruned exome), prepares Step 2 per-chromosome exome genotypes with updated sex, builds inputs (annotation/set/mask/covariates/phenotype), runs Step 1 and Step 2 (single-variant + gene-based burden/SKAT), then builds an HTML report and tables per phenotype.

#### **Config/Inputs**

+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| Key                                                | Description                                                                                                                                                                                | Required?                                        |
+====================================================+============================================================================================================================================================================================+==================================================+
| `project.name`                                     | Naming prefix included in all outputs                                                                                                                                                      | yes                                              |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `project.run_dir`                                  | Base directory for `preprocess/`, `input/`, `output/`, `logs/`                                                                                                                             | no (default `.`)                                 |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.vep_file`                                   | Combined VEP annotation CSV file (output from `select_variants_annotate.smk)`                                                                                                              | yes                                              |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.exome_chr`                                  | Per-chromosome exome pgen/psam path template (used for Step 2 genotypes, and Step 1 if `input.step1.input_type` is set to `exome`                                                          | yes                                              |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.pheno_file`                                 | Pre-built multi-phenotype TSV file creates separate outputs for each phenotype. Requires FID, IID phenotype columns. Defaults to single `STATUS` column built from cases/controls if empty | no                                               |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.phenotypes`                                 | Phenotype column names in `pheno_file` in report order                                                                                                                                     | no (default `[STATUS]`)                          |
|                                                    |                                                                                                                                                                                            |                                                  |
|                                                    | SIMPLEXO4: `[Overall, Under50, FamilyHistory, Malignant, ERPos, ERNeg]`                                                                                                                    |                                                  |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.step1.input_type`                           | Determines which type of data to use for Step 1.                                                                                                                                           | yes                                              |
|                                                    |                                                                                                                                                                                            |                                                  |
|                                                    | - `"array"` : uses PMBB's common SNPs file built from imputed genotype data (see PMBB's Freeze 4.0 Documentation)                                                                          |                                                  |
|                                                    |                                                                                                                                                                                            |                                                  |
|                                                    | - `"exome"` : uses PMBB's per-chromosome exome VCFs                                                                                                                                        |                                                  |
|                                                    |                                                                                                                                                                                            |                                                  |
|                                                    | SIMPLEXO4: used `array`                                                                                                                                                                    |                                                  |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.step1.array_all`                            | PLINK bfile prefix for the array QC step                                                                                                                                                   | yes if `input.step1.input_type` set to `array`   |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.cases`, `input.controls`                    | Sample list (one ID per line) for cases and controls, used to builds build the base sample set from these regardless                                                                       | yes                                              |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.step1_covariates`, `input.step2_covariates` | Covariate TSV files to be used in Step 1 and Step 2. Requires header; FID and IID as the first 2 columns.                                                                                  | yes                                              |
|                                                    |                                                                                                                                                                                            |                                                  |
|                                                    | SIMPLEXO4: `FID IID Age batch exome_PC1 exome_PC2 exome_PC3 exome_PC4 exome_PC5 exome_PC6`                                                                                                 |                                                  |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.sex_info`                                   | Sex TSV file without header to update genotypes. Three columns should be: FID, IID, sex                                                                                                    | yes                                              |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.mask_file`                                  | Custom mask definitions for regenie. First col is mask name, second is comma-separated components.                                                                                         | no                                               |
|                                                    |                                                                                                                                                                                            |                                                  |
|                                                    | SIMPLEXO4: Defaults to `M1 pathogenic` / `M2 pathogenic,vus`                                                                                                                               |                                                  |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.gene_consequence_exclude`                   | Specific exclusion rules dropped before pathogenic/VUS filtering. Formatted as GENE x consequence                                                                                          | no (default `mask_gene_consequence_exclude.txt`) |
|                                                    |                                                                                                                                                                                            |                                                  |
|                                                    | SIMPLEXO4: CHEK2 missense_variant                                                                                                                                                          |                                                  |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.negative_control_blacklist`                 | Variant IDs to scrub from the M4 synonymous negative-control mask                                                                                                                          | no                                               |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+
| `input.update_ids`                                 | ID remap file for `array_qc`, only used on the array Step 1 path. Formatted with first 2 cols as FID and IID of original, then next 2 cols as FID and IID to remap to.                     | no (default `input/array_update_ids.txt`)        |
+----------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+--------------------------------------------------+

#### **Submission**

``` bash
bsub -N -J array4 -eo logs/array4.e -oo logs/array4.o \
   snakemake -s scripts/regenie.smk --workflow-profile profiles/lsf -j 22 \
  --configfile config/array4.yaml
```

#### **Outputs**

- outputs/
  - `step2_single_variant_{PHENO}.regenie`, `step2_gene_based_{PHENO}.regenie` (per chromosome)
  - `chek2_carriers:` case/control counts for each variant in CHEK2
  - `consequence_matrix:` how many variants were in each Variant.Consequence category for each mask x thresh x gene
  - `gene_carriers:` case/control counts for each gene x mask x thresh for (top 50, suggestive ADD hits, FDR-significant ADD hits, fixed 15-gene list)
  - `examine_genes_ADD:` gene based test results for gene list to examine, with case/control counts
  - `top_genes_ADD` : gene based ADD test results for top 50 genes, with case/control counts
  - `top_genes_SKAT`: gene based SKAT test results for top 50 genes, with case/control counts
- reports/
  - `{PROJECT}_{PHENO}_report.html`: HTML summary report for each phenotype
  - `single_variant_results` : single-variant association analyses results
  - `significant_single_variants` : joined with VEP annotations + case/control counts
  - `gene_based_results_{ADDonly,alltests}` : gene-based burden test results
  - `CHEK2_M1_{ADDonly,alltests}` : gene-based burden results filtered down to M1, CHEK2 gene
  - `carriers_{topgenes,examinegenes}`: gene-based burden results with case/control counts for top 50 hits, then 15-gene examine list
  - `variant_contributions_{examinegenes,topgenes}` : individual variant contribution counts for top 50 hits, then 15-gene examine list

#### **Rules**

- `create_sample_list`
  - [input]{.underline}: cases, controls (config)
  - [output]{.underline}: {PROJECT}.samples.txt, {PROJECT}.samples_plink.txt
  - unions cases + controls into a single sample list, single-column (IID) and PLINK two-column (FID IID) formats
- `array_qc` (`array` Step 1 only)
  - [input]{.underline}: imputed array PLINK bfile (config), samples_plink.txt
  - [output]{.underline}: {PROJECT}.array.all_chr.step1.{pgen,pvar,psam,id,snplist}
  - QC's PMBB's imputed array data for use as Step 1 genotypes -- geno/MAF/HWE filters, sample ID remap (`--update-ids`)
- `get_exome_sample_ids`
  - [input]{.underline}: chr1 exome psam (config)
  - [output]{.underline}: {PROJECT}.exome.all_chr.step1.id
  - pulls the sample IDs present in the exome data so they can be intersected with the config/array sample lists below. Conditionally adds array IDs if step 1 is `array`
- `get_sample_list`
  - [input]{.underline}: exome sample IDs, config sample list, array sample IDs (if array-based Step 1)
  - [output]{.underline}: {PROJECT}.final_samples.txt
  - intersects exome sample IDs, config sample list, and conditionally adds array sample IDs into the final analysis sample set if step 1 is `array`
- `exome_chr_update_sex`
  - [input]{.underline}: per-chr exome pgen (config), final_samples.txt
  - [output]{.underline}: {PROJECT}.chr{CHR}.{pgen,pvar,psam}
  - updates sex from sex_info, restricts to the final sample set, applies a light per-chromosome genotype missingness filter (`--geno 0.1`) -- this is the Step 2 genotype set used directly by both step2 rules below
- `exome_chr_step1_filter` (`exome` Step 1 only)
  - [input]{.underline}: {PROJECT}.chr{CHR}.pgen (from exome_chr_update_sex)
  - [output]{.underline}: {PROJECT}.chr{CHR}.step1.{pgen,pvar,psam,snplist}
  - LD-prunes (`--indep-pairwise 1000 100 0.5`) and applies QC filters (MAF/geno/HWE), per chromosome
- `exome_merge_chr_step1` (exome Step 1 only)
  - [input]{.underline}: {PROJECT}.chr{CHR}.step1.pgen, all chromosomes
  - [output]{.underline}: {PROJECT}.exome.all_chr.step1.{pgen,pvar,psam}
  - merges all chromosomes' LD-pruned genotypes into one genome-wide Step 1 file
- `exome_merge_chr_snplist_step1` (exome Step 1 only)
  - [input]{.underline}: {PROJECT}.chr{CHR}.step1.snplist, all chromosomes
  - [output]{.underline}: {PROJECT}.exome.all_chr.step1.snplist
  - concatenates per-chromosome pruned SNP lists into one genome-wide list
- `preprocess_regenie`
  - [input]{.underline}: all_chr.vep.report.csv, samples.txt, cases/controls (or pre-built multi-phenotype file), step1/step2 covariate files
  - [output]{.underline}: {PROJECT}.regenie.{annotation,set,mask,covar.step1,covar.step2,pheno}.txt
  - calls `scripts/preprocess_regenie.py` to build all regenie-ready inputs (annotation/set/mask/covariate/phenotype files) from the VEP CSV
- `run_step1_regenie`
  - [input]{.underline}: Step 1 pgen/snplist (array or exome, per config), covar.step1, pheno
  - [output]{.underline}: {PROJECT}.step1_1.loco.gz, {PROJECT}.step1_pred.list
  - fits the regenie Step 1 null model
- `run_step2_single_variant`
  - [input]{.underline}: Step 2 per-chr pgen, covar.step2, pheno, step1_pred.list
  - [output]{.underline}: {PROJECT}.chr{CHR}.step2_single_variant\_{PHENO}.regenie
  - regenie Step 2 single-variant association test (Firth-corrected)
- `run_step2_gene_based`
  - [input]{.underline}: Step 2 per-chr pgen, covar.step2, pheno, annotation, set, mask, step1_pred.list
  - [output]{.underline}: {PROJECT}.chr{CHR}.step2_gene_based\_{PHENO}.regenie, {PROJECT}.chr{CHR}.step2_gene_based_masks.snplist
  - regenie Step 2 gene-based burden/SKAT test
    - `--write-mask-snplist` output (`masks_snplist`) is per-mask-definition, not per-phenotype, so it's one file per chromosome unlike every other output of this rule
- `create_vep_list`
  - [input]{.underline}: none (localrule)
  - [output]{.underline}: vep_files.list
  - wraps the single combined VEP CSV in the one-line-list format `mask_variant_stats.py` expects
- `mask_variant_stats`
  - [input]{.underline}: vep_files.list, pheno, per-chr pgen (all chromosomes)
  - [output]{.underline}: {PROJECT}.mask_variant_stats\_{PHENO}.tsv
  - per-mask, per-gene carrier/variant stats for the given phenotype, used to populate the report tables below
- `build_regenie_report_tables`
  - [input]{.underline}: step2_gene_based\_{PHENO}.regenie + masks.snplist (all chromosomes), mask_variant_stats\_{PHENO}.tsv, pheno
  - [output]{.underline}: {PROJECT}.report\_{PHENO}.{top_genes_ADD,top_genes_SKAT,variant_contrib_topgenes,variant_contrib_examinegenes,consequence_matrix,gene_carriers,chek2_carriers,examine_genes_ADD}.tsv
  - assembles the per-phenotype report tables used for inspection and to build the final report
- `generate_advanced_report`
  - [input]{.underline}: step2 single-variant + gene-based outputs (all chr), report tables, mask_variant_stats, pheno, report_regenie.Rmd
  - [output]{.underline}: reports/{PROJECT}\_{PHENO}\_report.html
  - renders the final per-phenotype HTML report

## Dependencies

+----------------------------------+-------------------------------------------------------------------+-----------------------------+
|                                  |                                                                   |                             |
+==================================+===================================================================+=============================+
| Script                           | **Purpose**                                                       | **Used by**                 |
+----------------------------------+-------------------------------------------------------------------+-----------------------------+
| `vep_vcf_parser.py`              | Parses VEP-annotated VCF -\> report CSV                           | select_variant_annotate.smk |
+----------------------------------+-------------------------------------------------------------------+-----------------------------+
| `preprocess_regenie.py`          | Builds annotation/set/mask/covariate/pheno files from the VEP CSV | regenie.smk                 |
+----------------------------------+-------------------------------------------------------------------+-----------------------------+
| `mask_variant_stats.py`          | Per-mask, per-gene carrier/variant stats for the report           | regenie.smk                 |
+----------------------------------+-------------------------------------------------------------------+-----------------------------+
| `build_regenie_report_tables.py` | Assembles top-genes/variant-contrib/carrier tables for the report | regenie.smk                 |
+----------------------------------+-------------------------------------------------------------------+-----------------------------+
| `report_regenie.Rmd`             | Renders the final per-phenotype HTML report                       | regenie.smk                 |
+----------------------------------+-------------------------------------------------------------------+-----------------------------+

#### vep_vcf_parser.py

- parses a VEP-annotated VCF into a flat per-variant CSV -- one row per variant, preferring a canonical + HIGH-impact annotation if one exists, else the first canonical transcript hit
- pulls each VEP core field and plugin/custom annotation block into columns: MANE, gnomAD per-population AFs, ClinVar/AutoGVP, SpliceAI, SIFT/PolyPhen/REVEL, AlphaMissense, MaveDB, LOFTEE
- computes `Variant.LoF_level` (the 1-4 pathogenicity tier) via `LofLevelAnnot` -- this is the column `preprocess_regenie.py` reads to build the pathogenic/VUS/synonymous mask categories
- `-m` mode controls sample-column handling: `no_sample` (drops all genotype/zygosity columns), `cohort`, `single,{sample.id}`, or `tumor_normal,{tumor.id},{normal.id}`
- optional `-g`/`-b`/`-r` gene-list/blacklist/BED-region filters, `-N` to drop caller columns, `-V` to include transformed VLR probability columns

#### preprocess_regenie.py

- takes the combined VEP CSV, sample/case/control lists (or a pre-built multi-phenotype file), and Step 1/Step 2 covariate files, and formats them for regenie's burden-test framework
- **annotation file**: one row per unique variant ID, most-pathogenic annotation kept on duplicates
  - `pathogenic`: Variant.LoF_level == 1
  - `vus`: Variant.LoF_level == 2
  - `synonymous`: negative-control category pulled from Variant.Consequence (M4), excludes anything already claimed by pathogenic/vus and anything in `--negative-control-blacklist`
- **set file**: for each gene, chr + position (from its first variant) + comma-separated variant IDs
- **mask file**: default is `M1 pathogenic` / `M2 pathogenic,vus`
  - Can be overridden with `--mask-file` (custom mask lines `Mx category1,category2,...`; categories are validated against {pathogenic, vus, synonymous} and checked for zero-variant categories)
- **covariate files**: passed from `--step1-covariates` / `--step2-covariates` (whatever columns they contain), filtered to the final sample set. If omitted, writes a minimal FID/IID-only file.
- **phenotype file**: single `STATUS` column (0/1) built from `--cases`/`--controls`, OR passed through from `--pheno-file` for multi-phenotype runs
- optional: `--gene-consequence-exclude` (drop specific GENE+consequence combos before pathogenic/VUS filtering)

#### mask_variant_stats.py

- reads every per-chromosome VEP report CSV (`--vep-list`) and keeps only mask-eligible sites. one row per unique variant ID, most-pathogenic kept on duplicates
- splits the phenotype file (`--pheno`/`--pheno-col`) into a cases-keep and a controls-keep file, then runs `plink2 --geno-counts` per chromosome against each, to get real per-variant carrier and called-sample counts (not REGENIE's fixed cohort-wide N)
- writes one TSV row per mask-eligible variant: `n_case_carrier`, `n_control_carrier`, `n_case_called`, `n_control_called`, plus derived `AC`/`AN`/`AF` -- this is the `--mask-stats` input that `build_regenie_report_tables.py` and `report_regenie.Rmd` both join against for real variant-level carrier counts

#### build_regenie_report_tables.py

- reads all per-chromosome `step2_gene_based_{PHENO}.regenie` + `masks.snplist` files and `mask_variant_stats.tsv`, and prebuilds the tables the Rmd renders from
- computes BH-adjusted FDR within test families (ADD alone; SKAT+SKATO together) and the top-N genes per mask×AAF-bin for both ADD and SKAT, and the consequence matrix
- `compute_gene_carriers()`: gets unique case/control carrier counts per gene×mask×bin through `plink2 --export A` dosage export (a sample counts if their dosage is \>0 for any variant in that gene's qualifying set) -- covers every key shown anywhere in the report: top-N, suggestive, FDR-significant, CHEK2 M1, and the fixed 15-gene list (`--examine-genes`)
- splits output into "top genes" vs "examine genes" (fixed watchlist) variant-contribution/carrier tables, since the Rmd consumes them separately

#### report_regenie.Rmd

- renders the final per-phenotype HTML report from the two Python scripts' prebuilt tables above, plus the `step2_single_variant`/`step2_gene_based` REGENIE files
- sections: project summary/QC (lambda GC, case/control counts), single-variant results (table + `Significant` flag), gene-based results (per-mask QQ/volcano plots, top-N tables, risk/protective FDR tables), directional burden, output file exports
- writes every `reports/` CSV/TSV export listed in the Outputs section above (`single_variant_results`, `significant_single_variants`, `gene_based_results_*`, `CHEK2_M1_*`, `carriers_*`, `variant_contributions_*`)
- render params (`project_name`, `data_dir`, `output_dir`, `pheno_file`, `pheno_col`, `step2_single_variant_mid`, `step2_gene_based_mid`) are passed in as absolute paths by `regenie.smk`'s `generate_advanced_report` rule

## 
