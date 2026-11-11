# From paired-end reads to sediment figures

Reproduce the sediment bacterial community analyses from paired FASTQs with two R scripts. `01_DADA2.R` writes `taxa.silva.xlsx`; after restarting R, `02_Sediments.R` reads that workbook, applies microDecon and produces the manuscript tables and figures. No precomputed phyloseq object, imputed spreadsheet or saved R workspace is an input. The authors ran both scripts on the study files and reproduced the reported numerical results and analysis figures.

## Files to obtain before running

Use both scripts and the input tables from the **same version** of [Project 2](https://github.com/CamiloMuMo/Bacterial-communities-of-Macrocystis-pyrifera). The [Releases page](https://github.com/CamiloMuMo/Bacterial-communities-of-Macrocystis-pyrifera/releases) is the distribution point for versioned code and control reads. The table below gives the file layout for the manuscript release.

 Input -> Where to get it and how it is used 
`01_DADA2.R`, `02_Sediments.R`, `samples.tsv`, `Env.csv` 
Place all four in the project root. `samples.tsv` is the sample-to-ENA map and controls the order of the libraries. `Env.csv` is the comma-separated table of **observed** metadata and geochemistry for the 31 environmental samples. Keep `NA` as missing; imputation is done only inside the secondary analysis.
Environmental reads: 25 sediment + 6 seawater pairs -> Use the 31 `ENA_run` accessions in `samples.tsv` to select runs from [ENA study PRJEB72006](https://www.ebi.ac.uk/ena/browser/view/PRJEB72006). Download both **raw, demultiplexed, primer-bearing** FASTQs for each run. Put `<ENA_run>_1.fastq.gz` in `R1/` and `<ENA_run>_2.fastq.gz` in `R2/`. Do not include kelp-tissue runs or already trimmed reads.
Control reads: `Blank1`, `Blank2`, `Blank3` (3 pairs) -> These controls have no ENA accessions. From the matching GitHub Release, obtain **both** `controls_fastq.zip` and `controls_fastq.zip.md5` as separate assets. The ZIP supplies the six original control FASTQs at the relative paths the script reads, plus a per-file checksum list. Verify them as shown below.
Taxonomy references -> Download the following two SILVA 138.2 files and place them in `Taxonomy/` without decompressing or renaming them:

- For `assignTaxonomy()` (classification through genus): [silva_nr99_v138.2_toGenus_trainset.fa.gz](https://zenodo.org/records/20955974/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1).
- For `addSpecies()` (species assignment): [silva_v138.2_assignSpecies.fa.gz](https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1). Do not substitute the `toSpecies_trainset` file for this step.

The checksums from the manuscript run are available in [reference_files.tsv](./Provenance%20Project%202/reference_files.tsv). Each new run of `01_DADA2.R` records the checksums of the files it used in `DADA2_output/reference_files.tsv`.
Leave the row order and the `ID` and `Library` columns in `samples.tsv` unchanged. The script accepts either the listed original FASTQ names (also with `.gz` appended) or the ENA filenames above for the 31 environmental runs. Use one copy of each mate. Controls are matched by their local filenames. Before running, there should be **34 libraries: 34 R1 files and 34 R2 files**, including the three control pairs. The six seawater samples are necessary to reconstruct the combined table and the descriptive phylum figure (Fig. S1); they are not replicates in the sediment tests.

## Verify and extract the controls

Start in a fresh project directory with the two release assets there. On Linux, run the following commands in **Bash**, not the R console:

```bash
md5sum -c controls_fastq.zip.md5 &&
  unzip -n controls_fastq.zip &&
  md5sum -c controls_MD5SUMS.txt
```

The first check verifies the downloaded ZIP **before** extraction. The last check must report `OK` for **six files** (three R1/R2 pairs). Stop if any check fails. The per-file MD5 values refer to the exact bytes distributed, including any existing `.gz` compression; do not decompress or recompress the FASTQs before checking. The checksum verifies the files, not the biological identity of a control.

## Software used

Install the packages listed at the top of each script before starting. Neither script installs software or changes library paths. The workflow used DADA2 **1.32.0**; the downstream analyses used R **4.5.3**, **microDecon 1.0.2** ([source commit `37307d8`](https://github.com/donaldtmcknight/microDecon/tree/37307d844384138fbc5a052e7f158a311b9f3291)), and missForest **1.6.1** with ranger **0.18.0**. Consult the two session reports for the R and package versions actually loaded in **each** run; use separate R environments if needed. No Excel application is required.

Step -> Session report written on each run 

DADA2 and workbook creation -> `DADA2_output/sessionInfo.txt`, `provenance/sessionInfo_01_DADA2.txt`
microDecon, ecological analyses and figures, `results/sessionInfo.txt`, `provenance/sessionInfo_02_Sediments.txt` 

These reports record the R version, operating system/platform, and versions of attached and namespace-loaded packages. They help diagnose differences between installations; they do not install packages, record the input files, or establish that the results match. Use the checksums, the ASV sequence map and the numerical checks below for those comparisons.

## Run in RStudio

Set the working directory to the extracted folder containing the two scripts (`Session > Set Working Directory > Choose Directory`). In a clean session run:

source("01_DADA2.R", encoding = "UTF-8")

The last step writes `taxa.silva.xlsx` directly into the working folder. The `counts` sheet has `OTU_ID`, then `Blank1`, `Blank2`, `Blank3`, then 31 environmental libraries. The `taxonomy` sheet contains the same ASV IDs, sequences and ranks. IDs are assigned in the order of the non-chimeric DADA2 sequence table, before sorting by nucleotide sequence to reproduce the historical merge. The exact sequence of each ASV (including `ASV_1`, `ASV_5` and `ASV_48`) is written to `DADA2_output/ASV_to_sequence.tsv` and `DADA2_output/ASV_sequences.fasta`. Copies supplied with the release are for comparison only. `DADA2_output/file_pairing.tsv` records which FASTQ files were actually used. An older one-sheet Excel file is not a substitute for this two-sheet workbook.

Restart R without restoring a saved workspace. From the same directory, run:

source("02_Sediments.R", encoding = "UTF-8")


Section 1 reads the counts worksheet and runs microDecon on the 34 libraries. The three blanks are controls for subtraction; they are not environmental replicates. Section 2 joins taxonomy and observed metadata by ID. Subsequent sections run alpha diversity, PERMANOVA/PERMDISP, NMDS, complete-case/partial/imputed RDA, prevalence-defined cores and taxonomic profiles. The 100-seed imputation check and 9,999-permutation tests can take time. The last sections export the figures and statistical tables. Existing outputs under these names are overwritten on rerun; use a fresh copy of the folder to keep a previous run.

Alternatively, run each command below from the project root in a terminal. Each starts a separate R session without restoring `.RData`:

bash
Rscript --vanilla 01_DADA2.R && Rscript --vanilla 02_Sediments.R


## Reading the outputs

Location -> Contents 

`results/supplement/S2A_*.tsv` -> Exact primer-position audit **before** trimming, pooled over the 31 environmental libraries; percentages use the total per direction.
`DADA2_output/after_trim_by_library.tsv` -> Separate diagnostic on the filtered reads: surviving shifted hits matched by read ID, observed residual prefixes, and full primers still identifiable in the first 60 nt. Ambiguous or mismatched primers are not ruled out. These data are not used to modify the denoising.
`results/supplement/S2B_*.tsv`, `S2C_*.tsv` -> Environmental merged lengths **before** chimera removal, and read tracking through the non-chimeric table; the controls are tracked separately. Detailed per-nucleotide lengths are in `DADA2_output/`.
`DADA2_output/error_learning_*.tsv` and `error_learning_*.txt` -> Eligible libraries, actual filtered bases and DADA2 error-learning output for this run. Both months, both materials and the three controls enter the joint workflow.
`results/supplement/` -> Manuscript tables: branch-specific ASV/read accounting (S2D–E), Hill numbers and exact tests (S3), community tests and permutation details (S4), RDA including site-by-date representation (S5), core (S6), focal sample (S7), and occurrence-filter RDA sensitivity (S8). Some panels contain more than one TSV.
`results/figures/` -> Editable SVGs: Figs. 2–6 and S1–S6 plus three additional sensitivity plots. BB_05 is labelled in the primary NMDS and RDA. See `results/tables/figure_manifest.tsv` for the exact filenames.
`results/tables/` -> Detailed outputs, including `Family_sample_profiles.tsv` (counts and percentages within each sediment), per-ASV scores, blank counts, model details and exact-test denominators.
`results/objects/` -> Newly computed phyloseq, analysis and figure objects, saved for later use; neither script requires these to exist.
`field/` -> The study's existing site map (Fig. 1), reported grain size (Table 2), and reported forest context. These are field or literature data, not generated from the sequence reads and not used as model inputs.

For statistical reporting, the S3–S5 exports give the sample sizes, degrees of freedom, statistics, R² where applicable, raw and adjusted P values and correction families for the relevant tests. Consult `results/supplement/S4_effective_permutations.tsv` for the effective permutation procedure in small groups, including the two August San Gregorio samples. The number of distinct group allocations is not necessarily the number of index or residual permutations. For a nonsignificant PERMDISP test, report “no heterogeneity detected”; equality of dispersions was not established.

Keep the two sensitivities distinct -> The primary bacterial matrix retains rare ASVs and all 25 sediments, including BB_05. In the occurrence sensitivity, detection in more than one of the **31** environmental libraries is evaluated before selecting sediments. Removing BB_05 is a separate change to the sample set. The three negative controls and the six seawater samples do not enter sediment community tests. The water data remain available for the descriptive phylum figure.

## Numbers for checking a local run

These are reference values from the revised analysis, **not** inputs read by the scripts. For Monte Carlo tests, sample order, package versions, permutation generation and seeds matter. An NMDS reflection or rotation by itself is not evidence that distances changed.

Check -> Reported value

Environmental reads at input; after filtering; merged; non-chimeric -> 5,915,385; 4,955,601; 4,526,601; 4,395,033 read pairs 
Pre-trim exact primer at nt 1 / shifted / no exact match (environmental reads, per direction) -> R1: 5,730,544 / 1,810 / 183,031; R2: 5,683,006 / 3,058 / 229,321 
Pre-chimera merged pairs of length 253 nt -> 4,281,360 of 4,526,601 (94.582%) 
Workbook before microDecon, including the controls | 24,899 ASVs × 34 libraries; non-chimeric blank reads: 20, 394 and 452 in Blank1–3 
Bacterial sediments, primary -> 25 samples; 20,429 ASVs; 3,361,111 reads 
Same 25 sediments with the occurrence filter -> 9,743 ASVs; 3,211,432 reads 
Primary matrix without BB_05 -> 24 samples; 20,363 ASVs; 3,235,552 reads 
Primary NMDS stress, with / without BB_05 -> 0.07150585 / 0.07462380 
March, common coverage, q = 1 (primary) -> Global exact P ≈ 0.0788; SG–BQ Holm P ≈ 0.0238; global gate closed 
Complete-case pooled RDA (16 samples) -> F = 2.464723, R² = 0.3812573, adjusted R² = 0.2265717, P = 0.0006 
Partial RDA, with / without BB_05 -Z Adjusted fraction ≈ 0.06283 / −0.00004; P = 0.0480 / 0.3922 
90% site core at detection ≥0.1% (primary) -> Bahía Buzo: 1; San Gregorio: 36; Buque Quemado: 80 

If raw, filtered or merged counts differ, check the library key, read pairing, raw-primer state and DADA2 settings first. If non-chimeric counts match but the bacterial table differs, compare ASV sequences, taxonomy reference hashes and control ordering before interpreting downstream P values. Do not change counts, ASV labels or thresholds to force agreement.

