# This file outlines how to run the PRS analysis. To perform quality control, imputation, and local ancestry analysis refer to: https://github.com/GP2code/SouthAfrican_PD_GWAS. The scripts included in this analysis were run on the HPC
# Pipeline overview is as follows:
1) MSP File Extraction & Covariate File Preparation:
	MSP files copied from GNOMIX output.
	extractHaploProportions.py used to infer ancestry.
	Age and sex merged from original GWAS covariates.
	Cleaned covariate file created (covar_prs.tsv).
2) Imputed VCF Processing:
	MAF filtered and converted to PLINK format.
	Merged into a single dataset, with related individuals removed.
3) Covariate File Matching:
	covar_prs_cleaned.txt matched to PLINK .fam IDs.
	Split into 70% (training/validation) and 30% (replication).
4) PRSice Setup & Execution:
	Loop to explore combinations of clump-kb, clump-r2, and clump-p.
	PRS run using Nalls and Kim summary stats on 70% subset.
	Optimal models evaluated using .summary files.
5) Testing Best Models:
	Best-performing parameters applied to 30% replication set.
6) AUC Evaluation:
	Phenotype recoded.
	AUC calculated via AUCBoot() and visualized with pROC.


######################################################################################
# Need covariate file. To do this use the MSP file from local ancestry analysis using gnmoix. Copy the msp files to the directory as follows
for dir in path_to_local_ancestry_output/STELLEN_*_NR; do
    if [[ -d "$dir" ]]; then  # Ensure it's a directory
        CHR=$(basename "$dir" | grep -oE '[0-9]+')  # Use -oE instead of -oP
        MSP_FILE="$dir/query_results.msp"

        if [[ -f "$MSP_FILE" ]]; then  # Check if the file exists
            cp "$MSP_FILE" "path_to_working_directory/msp_files/query_results_${CHR}.msp"
            echo "Copied $MSP_FILE to query_results_${CHR}.msp"
        else
            echo "Warning: \"$MSP_FILE\" not found in \"$dir\""
        fi
    fi
done

# To get the covar file with the inferred ancestry proportions run python script. We want a covariate file that has IID, AGE, SEX, ANCESTRY PROPORTIONS
# Run python script to process the fb files as follows:
python extractHaploProportions.py

# Add in age and sex to covar file using covar file from gwas
# Step 1: Strip header and sort both files by IID
tail -n +2 covarProjected_toModel.tsv | sort -k1,1 > covar_noheader_sorted.tsv
tail -n +2 proportionANCperIndividual.tsv | sort -k1,1 > props_noheader_sorted.tsv

awk -F'\t' '{OFS="\t"; print $1, $3, $4}' covar_noheader_sorted.tsv > covar_selected.tsv

# Sort both files by IID
sort -k1,1 covar_selected.tsv > covar_selected_sorted.tsv
sort -k1,1 proportionANCperIndividual.tsv > proportionANC_sorted.tsv

# Join the files by IID
join -t $'\t' -1 1 -2 1 covar_selected_sorted.tsv proportionANC_sorted.tsv > joined_tmp.tsv

# Add header to the final file
echo -e "IID\tSEX\tAGE\tAFR\tEUR\tMALAY\tNAMA\tSAS" > covar_prs.tsv
cat joined_tmp.tsv >> covar_prs.tsv


######################################################################################
# File prep for PRS
# Using the imputed files 
for chr in {1..22}; do
  # Filter by MAF >= 0.0005 and only biallelic SNPs
  bcftools view --min-af 0.0005 -Oz -o tmp_chr${chr}.biallelic.vcf.gz path_to_inputed_files/chr${chr}.dose.vcf.gz

  # Convert to plink format with missing variant IDs set
  plink2 --vcf tmp_chr${chr}.biallelic.vcf.gz dosage=DS \
    --set-all-var-ids @:#:\$r:\$a \
    --make-bed --out chr${chr}_plink
done

# Create a merge list
ls chr*_plink.bed | sed 's/.bed//' | grep -v "chr1_plink" > mergelist.txt

# Start with chr1_plink
plink --bfile chr1_plink --merge-list mergelist.txt --make-bed --out topmed_allchr_merged

# Remove related
plink --bfile topmed_allchr_merged --remove Related_individuals_toRemove.txt --make-bed --out all_chr_noRelated

# Rename so variant has the same IDs and column heading matches
# For Nalls sum stats (base file)
awk 'BEGIN {OFS="\t"} 
NR==1 {print $1, $2, $3, $4, "SNP", $5, $6, $7, $8, $9, $10; next} 
{print $1, $2, $3, $4, $1":"$2":"$4":"$3, $5, $6, $7, $8, $9, $10}' \
Nalls_fromGWASCatalog_GCST009325.tsv_hg38 > Nalls_hg38.txt

# Important! PRSice will only run when you have the same individuals as in your fam and covariate files. 
cut -d' ' -f2 all_chr_noRelated.fam | sort > fam_ids.txt
awk 'NR==FNR{a[$1]; next} $1 in a' fam_ids.txt covar_prs_cleaned.txt > covar_prs_matched.txt

# To run PRSice. Need to divide target cohort in 70% validation/hyperparameter optimisation and 30% testing.
cut -d' ' -f2 all_chr_noRelated.fam > all_ids.txt
shuf all_ids.txt > all_ids_shuffled.txt
total=$(wc -l < all_ids_shuffled.txt)
val_count=$(( (total * 70 + 99) / 100 ))  # Round up for 70%

head -n $val_count all_ids_shuffled.txt > ids_70_validation.txt
tail -n +$((val_count + 1)) all_ids_shuffled.txt > ids_30_test.txt

awk '{print 0, $1}' ids_70_validation.txt > ids_70_validation_fid0.txt
awk '{print 0, $1}' ids_30_test.txt > ids_30_test_fid0.txt

# For 70% validation (979 individuals; 445 cases; 534 controls)
plink --bfile all_chr_noRelated --keep ids_70_validation_fid0.txt --make-bed --out all_chr_70_validation

# For 30% test (419 individuals; 216 cases; 203 controls)
plink --bfile all_chr_noRelated --keep ids_30_test_fid0.txt --make-bed --out all_chr_30_test

# For 70% validation covariate file including the header
{
    printf "IID\tSEX\tAGE\tAFR\tEUR\tMALAY\tNAMA\tSAS\n"
    awk 'NR==FNR {a[$1]; next} $1 in a' ids_70_validation.txt covar_prs_cleaned.txt
} > covar_70_validation.txt

# For 30% test covariate file
{
    printf "IID\tSEX\tAGE\tAFR\tEUR\tMALAY\tNAMA\tSAS\n"
    awk 'NR==FNR {a[$1]; next} $1 in a' ids_30_test.txt covar_prs_cleaned.txt
} > covar_30_test.txt

# Prep summary stats for Kim
# Add SNP column to summary stats
awk 'BEGIN{OFS="\t"} NR==1{print $0, "SNP"} NR>1{print $0, $1 ":" $2 ":" $3 ":" $4}' Kim_2023_PD_MAMA.Joint.no23.hetFiltered.aggregate.trimmed.meta_hg38 > Kim_hg38.txt

sed '1s/(/_/g; 1s/)/_/g; 1s/BETA_FE/BETA/; 1s/\b\([A-Z_]*\)_\b/\1/g' Kim_hg38.txt > Kim_hg38_cleaned.txt


######################################################################################
# The following is for a status stratified approach (cases versus controls)

# Run on complete Nalls summary stats
# Initial run is as follows (https://github.com/choishingwan/PRSice), but you need to runs this for different values of --clump-kb, --clump-p, and --clump-r2. The prevalence flag also needs to be adjusted for phenotype prevalence. Run the validation cohort (70%) first to test parameters to get best modle then check how it performs in the independent replication cohort (30%). Discovery- summary stats: Base= summary stats
# For loop changing --clump-kb --clump-p and --clump-r2 (includes creating the output folder).

for kb in 100 250 500; do
  for r2 in 0.1 0.2 0.5 0.8; do
    for pval in 1e-3 1e-5 1e-6 5e-8; do
      outdir=results/LD_kb${kb}_r2${r2}_p${pval}
      mkdir -p ${outdir}
      outfile=${outdir}/NBA-PRS_kb${kb}_r2${r2}_p${pval}.summary

      if [[ -f "$outfile" ]]; then
        echo "Output exists for kb=${kb}, r2=${r2}, pval=${pval}. Skipping..."
        continue
      fi

      mkdir -p ${outdir}
      echo "Running PRSice with kb=${kb}, r2=${r2}, clump-p=${pval}"
      
      echo "Running PRSice with kb=${kb}, r2=${r2}, clump-p=${pval}"

      Rscript PRSice/PRSice.R --dir ./ \
        --prsice /apps/chpc/bio/PRSice/PRSice \
        --base summary_stats/Nalls_hg38.txt \
        --a1 effect_allele --a2 other_allele \
        --bp base_pair_location --snp SNP --chr chromosome \
        --pvalue p_value --beta \
        --target all_chr_70_validation --binary-target T \
        --clump-kb ${kb} \
        --clump-r2 ${r2} \
        --clump-p ${pval} \
        --interval 5e-6 --lower 0 --upper 0.5 \
        --num-auto 22 --seed 1234 --thread 1 --ignore-fid \
        --cov covar_70_validation.txt \
        --cov-col SEX,AGE,AFR,EUR,MALAY,NAMA,SAS \
        --perm 10000 --logit-perm --prevalence 0.1 \
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5 \
        --print-snp --quantile 10 --quant-ref 5 --quant-break 1,2,3,4,5,6,7,8,9,10 \
        --score sum \
        --out ${outdir}/NBA-PRS_kb${kb}_r2${r2}_p${pval}
    done
  done
done

# To determine which model ran the best use the .summary files in the results folder as the following:
echo -e "kb\tr2\tpval\tPRS_pval\tR2\tBeta" > prs_results_summary_1807.tsv

for kb in 100 250 500; do
  for r2 in 0.1 0.2 0.5 0.8; do
    for pval in 1e-3 1e-5 1e-6 5e-8; do
      file="results/LD_kb${kb}_r2${r2}_p${pval}/NBA-PRS_kb${kb}_r2${r2}_p${pval}.summary"
      if [[ -f $file ]]; then
        # Extract the header line number
        header_line=$(grep -n -m1 "^Phenotype" "$file" | cut -d: -f1)
        if [[ ! -z $header_line ]]; then
          # Get the first line after the header
          base_line=$(sed -n "$((header_line + 1))p" "$file")
          # Extract values from the Base line
          prs_r2=$(echo "$base_line" | cut -f4)
          prs_pval=$(echo "$base_line" | cut -f10)
          beta=$(echo "$base_line" | cut -f8)
          echo -e "${kb}\t${r2}\t${pval}\t${prs_pval}\t${prs_r2}\t${beta}" >> prs_results_summary_1807.tsv
        fi
      fi
    done
  done
done

# To find the best model run the following in R:
prs <- read.table("prs_results_summary_1807.tsv", header=TRUE, sep="\t")
prs$PRS_pval <- as.numeric(as.character(prs$PRS_pval))
prs$R2 <- as.numeric(as.character(prs$R2))
prs$Beta <- as.numeric(as.character(prs$Beta))

best_model <- prs[which.min(prs$PRS_pval), ]
print(best_model)

best_r2 <- prs[which.max(prs$R2), ]
print(best_r2)

# To plot this in R
library(ggplot2)
library(scales)  # for scientific notation in axis labels

# Load and prepare data
prs <- read.table("prs_results_summary.tsv", header=TRUE, sep="\t")
prs$PRS_pval <- as.numeric(as.character(prs$PRS_pval))
prs$R2 <- as.numeric(as.character(prs$R2))
prs$Beta <- as.numeric(as.character(prs$Beta))
prs$pval <- as.factor(prs$pval)  # for discrete x-axis

# 1. R² by p-value, colored by r2
ggplot(prs, aes(x=pval, y=R2, color=factor(r2), group=r2)) +
  geom_point(size=3) +
  geom_line() +
  facet_wrap(~kb, labeller=label_both) +
  scale_y_continuous(labels=scientific) +
  labs(
    title = "PRS R² across p-value thresholds",
    x = "Clumping p-value threshold",
    y = "R²",
    color = "r²"
  ) +
  theme_minimal(base_size = 14)

# 2. PRS p-value vs R²
ggplot(prs, aes(x=PRS_pval, y=R2, color=factor(kb), shape=factor(r2))) +
  geom_point(size=3) +
  scale_x_continuous(trans='log10', labels=scientific) +
  scale_y_continuous(labels=scientific) +
  labs(
    title = "PRS model trade-off: significance vs performance",
    x = "PRS p-value (log scale)",
    y = "R²",
    color = "kb",
    shape = "r²"
  ) +
  theme_minimal(base_size = 14)

# Based on the best model. Run the PRS again using the optimal inpute parameters determined above on the replication dataset (30%). 

for kb in 100; do
  for r2 in 0.5; do
    for pval in 1e-3; do
      outdir=results/30_dataset/LD_kb${kb}_r2${r2}_p${pval}
      mkdir -p ${outdir}
      
      echo "Running PRSice with kb=${kb}, r2=${r2}, clump-p=${pval}"

      Rscript PRSice/PRSice.R --dir ./ \
        --prsice /apps/chpc/bio/PRSice/PRSice \
        --base summary_stats/Nalls_hg38.txt \
        --a1 effect_allele --a2 other_allele \
        --bp base_pair_location --snp SNP --chr chromosome \
        --pvalue p_value --beta \
        --target all_chr_30_test --binary-target T \
        --clump-kb ${kb} \
        --clump-r2 ${r2} \
        --clump-p ${pval} \
        --interval 5e-6 --lower 0 --upper 0.5 \
        --num-auto 22 --seed 1234 --thread 1 --ignore-fid \
        --cov covar_30_test.txt \
        --cov-col SEX,AGE,AFR,EUR,MALAY,NAMA,SAS \
        --perm 10000 --logit-perm --prevalence 0.1 \
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5 \
        --print-snp --quantile 10 --quant-ref 5 --quant-break 1,2,3,4,5,6,7,8,9,10 \
        --score sum \
        --out ${outdir}/NBA-PRS_kb${kb}_r2${r2}_p${pval}
    done
  done
done


######################################################################################
# This needs to be repeated for the Kim summary statistics. Again start with the validation cohort (70%).

for kb in 100 250 500; do
  for r2 in 0.1 0.2 0.5 0.8; do
    for pval in 1e-3 1e-5 1e-6 5e-8; do
      outdir=results/Kim_prsice/70_dataset/LD_kb${kb}_r2${r2}_p${pval}
      mkdir -p ${outdir}
      
      outfile=${outdir}/NBA-PRS_kb${kb}_r2${r2}_p${pval}.summary

      if [[ -f "$outfile" ]]; then
        echo "Output exists for kb=${kb}, r2=${r2}, pval=${pval}. Skipping..."
        continue
      fi

      mkdir -p ${outdir}
      echo "Running PRSice with kb=${kb}, r2=${r2}, clump-p=${pval}"
      
     Rscript PRSice/PRSice.R --dir ./ \
        --prsice /apps/chpc/bio/PRSice/PRSice \
        --base summary_stats/Kim_hg38_cleaned.txt  \
        --a1 A1 --a2 A2 \
        --bp BP --snp SNP --chr CHR \
        --pvalue P_FE --beta \
        --target all_chr_70_validation --binary-target T \
        --clump-kb ${kb} \
        --clump-r2 ${r2} \
        --clump-p ${pval} \
        --interval 5e-6 --lower 0 --upper 0.5 \
        --num-auto 22 --seed 1234 --thread 1 --ignore-fid \
        --cov covar_70_validation.txt \
        --cov-col SEX,AGE,AFR,EUR,MALAY,NAMA,SAS \
        --perm 10000 --logit-perm --prevalence 0.1 \
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5 \
        --print-snp --quantile 10 --quant-ref 5 --quant-break 1,2,3,4,5,6,7,8,9,10 \
        --score sum \
        --out ${outdir}/NBA-PRS_kb${kb}_r2${r2}_p${pval}
    done
  done
done

# To determine which model ran the best use the .summary files in the Kim results folder as the following:
echo -e "kb\tr2\tpval\tPRS_pval\tR2\tBeta" > Kim_prs_results_summary.tsv

for kb in 100 250 500; do
  for r2 in 0.1 0.2 0.5 0.8; do
    for pval in 1e-3 1e-5 1e-6 5e-8; do
      file="results/Kim_prsice/70_dataset/LD_kb${kb}_r2${r2}_p${pval}/NBA-PRS_kb${kb}_r2${r2}_p${pval}.summary"
      if [[ -f $file ]]; then
        # Extract the header line number
        header_line=$(grep -n -m1 "^Phenotype" "$file" | cut -d: -f1)
        if [[ ! -z $header_line ]]; then
          # Get the first line after the header
          base_line=$(sed -n "$((header_line + 1))p" "$file")
          # Extract values from the Base line
          prs_r2=$(echo "$base_line" | cut -f4)
          prs_pval=$(echo "$base_line" | cut -f10)
          beta=$(echo "$base_line" | cut -f8)
          echo -e "${kb}\t${r2}\t${pval}\t${prs_pval}\t${prs_r2}\t${beta}" >> Kim_prs_results_summary.tsv
        fi
      fi
    done
  done
done

# To find the best model run the following in R:
prs <- read.table("Kim_prs_results_summary.tsv", header=TRUE, sep="\t")
prs$PRS_pval <- as.numeric(as.character(prs$PRS_pval))
prs$R2 <- as.numeric(as.character(prs$R2))
prs$Beta <- as.numeric(as.character(prs$Beta))

best_model <- prs[which.min(prs$PRS_pval), ]
print(best_model)

best_r2 <- prs[which.max(prs$R2), ]
print(best_r2)
      
# Based on the best model. Run the PRS again using the optimal inpute parameters determined above on the replication dataset (30%).     

for kb in 100; do
  for r2 in 0.5; do
    for pval in 1e-3; do
      outdir=results/Kim_prsice/30_dataset/LD_kb${kb}_r2${r2}_p${pval}
      mkdir -p ${outdir}
      echo "Running PRSice with kb=${kb}, r2=${r2}, clump-p=${pval}"

     Rscript PRSice/PRSice.R --dir ./ \
        --prsice /apps/chpc/bio/PRSice/PRSice \
        --base summary_stats/Kim_hg38_cleaned.txt  \
        --a1 A1 --a2 A2 \
        --bp BP --snp SNP --chr CHR \
        --pvalue P_FE --beta \
        --target all_chr_30_test --binary-target T \
        --clump-kb ${kb} \
        --clump-r2 ${r2} \
        --clump-p ${pval} \
        --interval 5e-6 --lower 0 --upper 0.5 \
        --num-auto 22 --seed 1234 --thread 1 --ignore-fid \
        --cov covar_30_test.txt \
        --cov-col SEX,AGE,AFR,EUR,MALAY,NAMA,SAS \
        --perm 10000 --logit-perm --prevalence 0.1 \
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5 \
        --print-snp --quantile 10 --quant-ref 5 --quant-break 1,2,3,4,5,6,7,8,9,10 \
        --score sum \
        --out ${outdir}/NBA-PRS_kb${kb}_r2${r2}_p${pval}
    done
  done
done


######################################################################################
# PRSice2 does not calculate AUC, which is needed to make PRS approaches comparable. To calculate it run the following using the .best output file from the 30 split from PRSice. To understand how the AUCBoot script works you can refer to AUCBoot_codes.R. Need to do done for both Nalls and Kim summary stats. 

# Run in R
library(bigstatsr)

# Load fam file, no header
fam_df <- read.table("all_chr_30_test.fam", header = FALSE)

# Assign column names
colnames(fam_df) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")

# Check unique phenotype values
table(fam_df$PHENO)

# Load your PRS .best file (adjust the path/filename!)
prs_df <- read.table("NBA-PRS_kb100_r20.5_p1e-3.best", header = TRUE)

# Merge with PRS file by FID and IID
merged_df <- merge(prs_df, fam_df[, c("FID", "IID", "PHENO")], by = c("FID", "IID"))

# Recode phenotype: 2 = case -> 1, 1 = control -> 0
phenotype <- ifelse(merged_df$PHENO == 2, 1, 0)

# PRS values
prs_values <- merged_df$PRS

# Check lengths to make sure they match
length(phenotype)
length(prs_values)

# Run AUCBoot
AUCBoot(pred = prs_values, target = phenotype, nboot = 10000, seed = 123)

# To plot AUC 
# Run in R
install.packages("pROC")  # run this if you haven't installed pROC
library(pROC)

# Create an ROC curve object:
roc_obj <- roc(response = phenotype, predictor = prs_values)

# Plot the ROC curve:
plot(roc_obj, col = "#098dc3ff", lwd = 2, main = "ROC Curve for PRS Model using Kim et al 2024")

# Add AUC to the plot (optional):
auc_val <- auc(roc_obj)
legend("bottomright", legend = paste("AUC =", round(auc_val, 3)), col = "#098dc3ff", lwd = 2)

# Set up PNG device with 300 DPI and larger size for quality
png("KIM_PRS_ROC_Curve_300DPI.png", width = 6*300, height = 6*300, res = 300)  # 6x6 inches at 300 dpi

# Plot ROC curve with specified color and title
plot(roc_obj, col = "#33a02c", lwd = 2, main = "ROC Curve for PRS model using Kim et al 2023")

# Add AUC legend
legend("bottomright", legend = paste("AUC =", round(auc(roc_obj), 3)), col = "#33a02c", lwd = 2)

# Close the device
dev.off()


### To plot both AUC
# Load necessary library
library(pROC)

# --- File paths ---
fam_path <- "all_chr_30_test.fam"
nalls_best <- "nalls_NBA-PRS_kb100_r20.5_p1e-3.best"
kim_best <- "kim_NBA-PRS_kb100_r20.5_p1e-3.best"


# --- Load .fam file ---
fam_df <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
colnames(fam_df) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
fam_df$PHENO <- ifelse(fam_df$PHENO == 2, 1, 0)  # 2=case, 1=control → recode to 1=case, 0=control

# --- Load .best files ---
nalls_df <- read.table(nalls_best, header = TRUE, stringsAsFactors = FALSE)
kim_df <- read.table(kim_best, header = TRUE, stringsAsFactors = FALSE)

# --- Merge PRS with phenotype ---
nalls_merged <- merge(nalls_df, fam_df[, c("FID", "IID", "PHENO")], by = c("FID", "IID"))
kim_merged <- merge(kim_df, fam_df[, c("FID", "IID", "PHENO")], by = c("FID", "IID"))

# --- Prepare vectors ---
prs_nalls <- nalls_merged$PRS
pheno_nalls <- nalls_merged$PHENO

prs_kim <- kim_merged$PRS
pheno_kim <- kim_merged$PHENO

# --- Compute ROC curves ---
roc_nalls <- roc(pheno_nalls, prs_nalls, quiet = TRUE)
roc_kim <- roc(pheno_kim, prs_kim, quiet = TRUE)

# --- Save high-res ROC plot ---
png(filename = "PRS_ROC_Comparison.png", width = 6, height = 6, units = "in", res = 400)
plot(roc_nalls, col = "#1f78b4", lwd = 3,
     main = "ROC Curve: PRS Models (Nalls et al. vs Kim et al.)",
     legacy.axes = TRUE, cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.1)
plot(roc_kim, col = "#33a02c", lwd = 3, add = TRUE)

legend("bottomright",
       legend = c(sprintf("Nalls et al. 2019 (AUC = %.3f)", auc(roc_nalls)),
                  sprintf("Kim et al. 2024 (AUC = %.3f)", auc(roc_kim))),
       col = c("#1f78b4", "#33a02c"),
       lwd = 3,
       cex = 1.1)

dev.off()

# --- Optional: Print AUC values to console ---
cat(sprintf("AUC Nalls: %.3f\n", auc(roc_nalls)))
cat(sprintf("AUC Kim: %.3f\n", auc(roc_kim)))


######################################################################################
# To look at the SNPs contributing most to the predictive model. You will need to run PRSice with --print-snp and --score flags to be able to do this section. Run the top_snps_prsice.R script on the highest performing PRSice input for the 70%. You'll need to run it for Nalls and for Kim. 

# For Nalls
# From the .snp output file, filter the variants that match your p-value. this output file should have the same number of snps as indicated in your .summary file from your PRSice2 run. The threshold/ cut off used in the filtering is based on the value in the Threshold column of the .summary file 
awk '$4 < 0.000075' NBA-PRS_kb100_r20.5_p1e-3.snp > nalls_snps.txt

# Then need to run PRSice2 looped removing one SNP at a time:
# Step 1: Create a filtered base file with only the SNPs in the SNP_LIST
FILTERED_BASE="tmp/Nalls_hg38_filtered.txt"
SNP_LIST="nalls_snps.txt"
BASE_FILE="summary_stats/Nalls_hg38.txt"
OUTDIR="variants/nalls"
mkdir -p tmp "${OUTDIR}"

awk 'NR==FNR {snps[$2]; next} FNR==1 || ($5 in snps)' ${SNP_LIST} ${BASE_FILE} > ${FILTERED_BASE}

# Step 2: Loop through each SNP in the filtered base file (column 5)
awk 'NR > 1 {print $5}' ${FILTERED_BASE} | while read -r snp; do
echo "Running PRSice with kb=100, r2=0.5, clump-p=1e-3, excluding SNP=${snp}"

 # Create a temporary base file excluding the current SNP
# tmp_base="tmp/Nalls_hg38_filtered_excl_${snp//:/_}.txt"
awk -v exclude="$snp" 'BEGIN {FS=OFS="\t"} FNR==1 || $5 != exclude' ${FILTERED_BASE} > tmp/Nalls_hg38_filtered_excl_${snp//:/_}.txt

        # Run PRSice
        Rscript PRSice/PRSice.R --dir ./ \
          --prsice /apps/chpc/bio/PRSice/PRSice \
          --base ${tmp_base} \
          --a1 effect_allele --a2 other_allele \
          --bp base_pair_location --snp SNP --chr chromosome \
          --pvalue p_value --beta \
          --target all_chr_70_validation --binary-target T \
          --clump-kb 100 \
          --clump-r2 0.5 \
          --clump-p 1e-3 \
          --interval 5e-6 --lower 0 --upper 0.5 \
          --num-auto 22 --seed 1234 --thread 1 --ignore-fid \
          --cov covar_70_validation.txt \
          --cov-col SEX,AGE,AFR,EUR,MALAY,NAMA,SAS \
          --perm 10000 --logit-perm --prevalence 0.1 \
          --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5 \
          --print-snp --quantile 10 --quant-ref 5 --quant-break 1,2,3,4,5,6,7,8,9,10 \
          --score sum \
          --out ${outdir}/PRS_kb100_r20.5_p1e-3_noSNP_${snp//:/_}

        # Optionally delete temp base file
#        rm ${tmp_base}

 done

# The variant-containing window was extracted from the imputed files and matched it to the inferred ancestry in the marginal segment probability (MSP) file. Refer to the https://github.com/GP2code/SouthAfrican_PD_GWAS scripts to see how to locate the local ancestry window

# Run the following in R
library(pROC)

# --- File paths ---
fam_path <- "all_chr_70_validation.fam"
best_dir <- "variants/nalls"
output_csv <- "Nalls_AUC_per_excluded_snp.csv"
nalls_top_contributors_csv <- "Nalls_Top_SNP_Contributors.csv"
nalls_top_predictive_csv <- "Nalls_Top_10_Predictive_SNPs.csv"

# --- Original full-model AUC --- get this from where you originally calculated AUC
original_auc <- 0.6287

# --- Load .fam file ---
fam_df <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
colnames(fam_df) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
fam_df$PHENO <- ifelse(fam_df$PHENO == 2, 1, 0)  # recode: 2=case → 1, 1=control → 0

# --- Get all .best files ---
best_files <- list.files(path = best_dir, pattern = "*.best$", full.names = TRUE)

# --- Create results container ---
results <- data.frame(SNP = character(), AUC = numeric(), stringsAsFactors = FALSE)

# --- Loop through files and calculate AUC ---
for (file in best_files) {
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  merged <- merge(df, fam_df[, c("FID", "IID", "PHENO")], by = c("FID", "IID"))

  if (nrow(merged) > 0) {
    roc_obj <- roc(merged$PHENO, merged$PRS, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))

    # Extract SNP name from filename
    fname <- basename(file)
    snp <- sub(".*noSNP_(.*)\\.best", "\\1", fname)
    snp <- gsub("_", ":", snp)

    results <- rbind(results, data.frame(SNP = snp, AUC = auc_val))
  }
}

# --- Calculate AUC difference from original model ---
results$delta_AUC <- results$AUC - original_auc

# --- Save full results ---
write.csv(results, output_csv, row.names = FALSE)
cat("All AUC results saved to", output_csv, "\n")

# --- Top 10 SNPs that had the greatest overall impact (absolute delta) ---
top_contributors <- results[order(abs(results$delta_AUC), decreasing = TRUE), ]
write.csv(head(top_contributors, 10), nalls_top_contributors_csv, row.names = FALSE)
cat("Top 10 highest-impact SNPs saved to", nalls_top_contributors_csv, "\n")

# --- Top 10 SNPs that contributed most positively (most drop in AUC when removed) ---
top_predictive <- results[order(results$delta_AUC), ]  # most negative delta first
write.csv(head(top_predictive, 10), nalls_top_predictive_csv, row.names = FALSE)
cat("Top 10 predictive SNPs saved to", nalls_top_predictive_csv, "\n")

# For Kim
# From the .snp output file, filter the variants that match your p-value. this output file should have the same number of snps as indicated in your .summary file from your PRSice2 run. The threshold/ cut off used in the filtering is based on the value in the Threshold column of the .summary file 
awk '$4 < 0.000345' NBA-PRS_kb100_r20.5_p1e-3.snp > kim_snps.txt

# Then need to run PRSice2 looped removing one SNP at a time:
# Paths
BASE_FILE="summary_stats/Kim_hg38_cleaned.txt"
SNP_LIST="kim_snps.txt"
FILTERED_BASE="tmp/kim/Kim_hg38_filtered.txt"
OUTDIR="variants/kim"

# Step 2: Loop through each SNP in the filtered base file (column 5 = SNP)
awk 'NR > 1 {print $14}' ${FILTERED_BASE} | while read -r snp; do
    clean_snp=${snp//:/_}
    output_summary="${OUTDIR}/PRS_kb100_r20.5_p1e-3_noSNP_${clean_snp}.summary"

    if [[ -f "${output_summary}" ]]; then
        echo "Skipping ${snp}, output already exists."
        continue
    fi

    echo "Running PRSice excluding SNP=${snp}, kb=100, r2=0.5, pval=1e-3"

    # Create temporary base file excluding the SNP
    tmp_base="tmp/kim/Kim_hg38_filtered_excl_${clean_snp}.txt"
    awk -v exclude="$snp" 'BEGIN {FS=OFS="\t"} FNR==1 || $14 != exclude' ${FILTERED_BASE} > "${tmp_base}"

    # Run PRSice
    Rscript PRSice/PRSice.R --dir ./ \
      --prsice /apps/chpc/bio/PRSice/PRSice \
      --base "${tmp_base}" \
      --a1 A1 --a2 A2 \
      --bp BP --snp SNP --chr CHR \
      --pvalue P_FE --beta \
      --target all_chr_70_validation --binary-target T \
      --clump-kb 100 \
      --clump-r2 0.5 \
      --clump-p 1e-3 \
      --interval 5e-6 --lower 0 --upper 0.5 \
      --num-auto 22 --seed 1234 --thread 1 --ignore-fid \
      --cov covar_70_validation.txt \
      --cov-col SEX,AGE,AFR,EUR,MALAY,NAMA,SAS \
      --perm 10000 --logit-perm --prevalence 0.1 \
      --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5 \
      --print-snp --quantile 10 --quant-ref 5 --quant-break 1,2,3,4,5,6,7,8,9,10 \
      --score sum \
      --out "${OUTDIR}/PRS_kb100_r20.5_p1e-3_noSNP_${clean_snp}"

    # Clean up (optional)
    # rm -f "${tmp_base}"

done

# Once this is completed, run the following in R to determine the variants contributing most to the model. You can then cross reference the variant carries with their local ancestry information to confirm the local ancestry window. 
library(pROC)

# --- File paths ---
fam_path <- "all_chr_70_validation.fam"
best_dir <- "variants/kim"
output_csv <- "Kim_AUC_per_excluded_snp.csv"
nalls_top_contributors_csv <- "Kim_Top_SNP_Contributors.csv"
nalls_top_predictive_csv <- "Kim_Top_10_Predictive_SNPs.csv"

# --- Original full-model AUC ---
original_auc <- 0.6294

# --- Load .fam file ---
fam_df <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE)
colnames(fam_df) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
fam_df$PHENO <- ifelse(fam_df$PHENO == 2, 1, 0)  # recode: 2=case → 1, 1=control → 0

# --- Get all .best files ---
best_files <- list.files(path = best_dir, pattern = "*.best$", full.names = TRUE)

# --- Create results container ---
results <- data.frame(SNP = character(), AUC = numeric(), stringsAsFactors = FALSE)

# --- Loop through files and calculate AUC ---
for (file in best_files) {
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  merged <- merge(df, fam_df[, c("FID", "IID", "PHENO")], by = c("FID", "IID"))

  if (nrow(merged) > 0) {
    roc_obj <- roc(merged$PHENO, merged$PRS, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))

    # Extract SNP name from filename
    fname <- basename(file)
    snp <- sub(".*noSNP_(.*)\\.best", "\\1", fname)
    snp <- gsub("_", ":", snp)

    results <- rbind(results, data.frame(SNP = snp, AUC = auc_val))
  }
}

# --- Calculate AUC difference from original model ---
results$delta_AUC <- results$AUC - original_auc

# --- Save full results ---
write.csv(results, output_csv, row.names = FALSE)
cat("All AUC results saved to", output_csv, "\n")

# --- Top 10 SNPs that had the greatest overall impact (absolute delta) ---
top_contributors <- results[order(abs(results$delta_AUC), decreasing = TRUE), ]
write.csv(head(top_contributors, 10), kim_top_contributors_csv, row.names = FALSE)
cat("Top 10 highest-impact SNPs saved to", kim_top_contributors_csv, "\n")

# --- Top 10 SNPs that contributed most positively (most drop in AUC when removed) ---
top_predictive <- results[order(results$delta_AUC), ]  # most negative delta first
write.csv(head(top_predictive, 10), kim_top_predictive_csv, row.names = FALSE)
cat("Top 10 predictive SNPs saved to", kim_top_predictive_csv, "\n")


#####################################################################################
# Run then to see variance by covariates: For each run remove/ change the combination of covariates included in the analysis after the --cov-col flag. You want to run the following seven combinations: 1. age, 2. sex, 3. ancestries, 4. age, sex, 5. age, ancestries, 6. sex, ancestries, 7. age, sex, ancestries

# The run here will be using the Nalls summary statistics with the validation cohort (70%)
for kb in 100; do
  for r2 in 0.5; do
    for pval in 1e-3; do
      outdir=covar_results/nalls/70_dataset/sex
      mkdir -p ${outdir}
      
      echo "Running PRSice with kb=${kb}, r2=${r2}, clump-p=${pval}"

      Rscript PRSice/PRSice.R --dir ./ \
        --prsice /apps/chpc/bio/PRSice/PRSice \
        --base summary_stats/Nalls_hg38.txt \
        --a1 effect_allele --a2 other_allele \
        --bp base_pair_location --snp SNP --chr chromosome \
        --pvalue p_value --beta \
        --target all_chr_70_validation --binary-target T \
        --clump-kb ${kb} \
        --clump-r2 ${r2} \
        --clump-p ${pval} \
        --interval 5e-6 --lower 0 --upper 0.5 \
        --num-auto 22 --seed 1234 --thread 1 --ignore-fid \
        --cov covar_70_validation.txt \
        --cov-col SEX \
        --perm 10000 --logit-perm --prevalence 0.1 \
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5 \
        --print-snp --quantile 10 --quant-ref 5 --quant-break 1,2,3,4,5,6,7,8,9,10 \
        --score sum \
        --out ${outdir}
    done
  done
done

# The next run here will be using the Nalls summary statistics with the replication cohort (30%)
for kb in 100; do
  for r2 in 0.5; do
    for pval in 1e-3; do
      outdir=covar_results/nalls/30_dataset/sex
      mkdir -p ${outdir}
      
      echo "Running PRSice with kb=${kb}, r2=${r2}, clump-p=${pval}"

      Rscript PRSice/PRSice.R --dir ./ \
        --prsice /apps/chpc/bio/PRSice/PRSice \
        --base summary_stats/Nalls_hg38.txt \
        --a1 effect_allele --a2 other_allele \
        --bp base_pair_location --snp SNP --chr chromosome \
        --pvalue p_value --beta \
        --target all_chr_30_test --binary-target T \
        --clump-kb ${kb} \
        --clump-r2 ${r2} \
        --clump-p ${pval} \
        --interval 5e-6 --lower 0 --upper 0.5 \
        --num-auto 22 --seed 1234 --thread 1 --ignore-fid \
        --cov covar_30_test.txt \
        --cov-col SEX \
        --perm 10000 --logit-perm --prevalence 0.1 \
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5 \
        --print-snp --quantile 10 --quant-ref 5 --quant-break 1,2,3,4,5,6,7,8,9,10 \
        --score sum \
        --out ${outdir}
    done
  done
done

# The next run here will be using the Kim summary statistics with the validation cohort (70%):
for kb in 100 250 500; do
  for r2 in 0.1 0.2 0.5 0.8; do
    for pval in 1e-3 1e-5 1e-6 1e-8; do
      outdir=results/Kim_prsice/70_dataset/LD_kb${kb}_r2${r2}_p${pval}
      mkdir -p ${outdir}
      
      echo "Running PRSice with kb=${kb}, r2=${r2}, clump-p=${pval}"

     Rscript PRSice/PRSice.R --dir ./ \
        --prsice /apps/chpc/bio/PRSice/PRSice \
        --base summary_stats/Kim_hg38_cleaned.txt  \
        --a1 A1 --a2 A2 \
        --bp BP --snp SNP --chr CHR \
        --pvalue P_FE --beta \
        --target all_chr_70_validation --binary-target T \
        --clump-kb ${kb} \
        --clump-r2 ${r2} \
        --clump-p ${pval} \
        --interval 5e-6 --lower 0 --upper 0.5 \
        --num-auto 22 --seed 1234 --thread 1 --ignore-fid \
        --cov covar_70_validation.txt \
        --cov-col SEX,AGE,AFR,EUR,MALAY,NAMA,SAS \
        --perm 10000 --logit-perm --prevalence 0.1 \
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5 \
        --print-snp --quantile 10 --quant-ref 5 --quant-break 1,2,3,4,5,6,7,8,9,10 \
        --score sum \
        --out ${outdir}/NBA-PRS_kb${kb}_r2${r2}_p${pval}
    done
  done
done

# The next run here will be using the Kim summary statistics with the replication cohort (30%):
for kb in 100; do
  for r2 in 0.5; do
    for pval in 1e-3; do
      outdir=results/Kim_prsice/30_dataset/LD_kb${kb}_r2${r2}_p${pval}
      mkdir -p ${outdir}
      echo "Running PRSice with kb=${kb}, r2=${r2}, clump-p=${pval}"

     Rscript PRSice/PRSice.R --dir ./ \
        --prsice /apps/chpc/bio/PRSice/PRSice \
        --base summary_stats/Kim_hg38_cleaned.txt  \
        --a1 A1 --a2 A2 \
        --bp BP --snp SNP --chr CHR \
        --pvalue P_FE --beta \
        --target all_chr_30_test --binary-target T \
        --clump-kb ${kb} \
        --clump-r2 ${r2} \
        --clump-p ${pval} \
        --interval 5e-6 --lower 0 --upper 0.5 \
        --num-auto 22 --seed 1234 --thread 1 --ignore-fid \
        --cov covar_30_test.txt \
        --cov-col SEX,AGE,AFR,EUR,MALAY,NAMA,SAS \
        --perm 10000 --logit-perm --prevalence 0.1 \
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5 \
        --print-snp --quantile 10 --quant-ref 5 --quant-break 1,2,3,4,5,6,7,8,9,10 \
        --score sum \
        --out ${outdir}/NBA-PRS_kb${kb}_r2${r2}_p${pval}
    done
  done
done


######################################################################################
# For additional sensitivity analysis run the following on the .best output files in R: 
library(data.table)
library(pROC)
library(caret)

# Input paths
prs_file <- "nalls_70_NBA-PRS_kb100_r20.5_p1e-3.best"
fam_file <- "all_chr_70_validation.fam"

# Load PRS
prs <- fread(prs_file)

# Load .fam file (no header)
fam <- fread(fam_file, header = FALSE)
setnames(fam, c("FID", "IID", "PID", "MID", "SEX", "PHENO"))

# Recode phenotype (2 = case → 1, 1 = control → 0, else NA)
fam[, CASE := fifelse(PHENO == 2, 1, fifelse(PHENO == 1, 0, NA_integer_))]

# Remove individuals with missing phenotype
fam <- fam[!is.na(CASE)]

# Merge PRS with phenotype
dat <- merge(prs, fam[, .(FID, IID, CASE)], by = c("FID", "IID"))

# Standardize PRS using control mean and SD
mean_controls <- mean(dat$PRS[dat$CASE == 0])
sd_controls <- sd(dat$PRS[dat$CASE == 0])
dat[, zSCORE := (PRS - mean_controls) / sd_controls]

# Logistic regression
model <- glm(CASE ~ zSCORE, family = "binomial", data = dat)
dat[, probDisease := predict(model, dat, type = "response")]

# ROC
dat[, reported := ifelse(CASE == 1, "DISEASE", "CONTROL")]
roc_obj <- roc(response = dat$reported, predictor = dat$probDisease, quiet = TRUE)
auc_val <- auc(roc_obj)

# Fixed threshold
threshold <- 0.5

# Predicted classes
dat[, predicted := ifelse(probDisease > threshold, "DISEASE", "CONTROL")]

# Confusion matrix
conf_mat <- confusionMatrix(
  factor(dat$predicted, levels = c("CONTROL", "DISEASE")),
  factor(dat$reported, levels = c("CONTROL", "DISEASE")),
  positive = "DISEASE"
)

# Output results
results_dt <- data.table(
  Model = "nalls_70",
  AUC = auc_val,
  Accuracy = conf_mat$overall["Accuracy"],
  CI_Lower = conf_mat$overall["AccuracyLower"],
  CI_Upper = conf_mat$overall["AccuracyUpper"],
  Balanced_Accuracy = conf_mat$byClass["Balanced Accuracy"],
  Sensitivity = conf_mat$byClass["Sensitivity"],
  Specificity = conf_mat$byClass["Specificity"]
)

# Save results
fwrite(results_dt, "nalls_70_model2_results.tsv", sep = "\t")
print(results_dt)
