suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(sva)))

args <- commandArgs(trailingOnly = TRUE)
# $1 = plink_pref
# $2 = metadata/covariates file (csv)
# $3 = workdir
# $4 = directory the script was run from (for gene_loc.txt)
# $5 = single cell data name (for expression matrices)

workdir <- args[[3]]
plink_pref <- args[[1]]


gene_annotation <- read.table(paste0(args[[4]], "/gene_loc.txt"), header = TRUE) %>%
  distinct(NAME, .keep_all = TRUE)

wgs_subset <- read.table(paste0(plink_pref, ".fam"), header = FALSE)

genotype_pcs <- read.table(paste0(workdir, "/", basename(plink_pref), ".eigenvec")) %>%
  rename(
    Sample = V2,
    G_PC1 = V3,
    G_PC2 = V4,
    G_PC3 = V5,
    G_PC4 = V6,
    G_PC5 = V7
  ) %>%
  select(!V1)

metadata <- read.csv(args[[2]]) %>%
  mutate(Sex=as.factor(Sex), Age=as.numeric(gsub("[^0-9.-]", "", Age))) %>%
  select(Sample, Sex, Age)


file <- args[[5]]

# Subset only the NPH expression and composition files
efile <- paste0(workdir, "/processed_matrices/", file, "_expression_matrix_ds.csv")
cfile <- paste0(workdir, "/processed_matrices/", file, "_composition_matrix_ds.csv")

# Read expression data
edata <- t(
  as.matrix(read.csv(efile) %>% 
    column_to_rownames("X"))
)
colnames(edata) <- sub("_NPH$", "", colnames(edata))
# Filter metadata + expression data to match genotype data
pheno <- metadata %>%
  filter(Sample %in% wgs_subset$V2 & Sample %in% colnames(edata))  %>%
  arrange(match(Sample, wgs_subset$V2))
rownames(pheno) <- pheno$Sample
edata <- edata[, colnames(edata) %in% wgs_subset$V2 & colnames(edata) %in% rownames(pheno)]

# 1. Align Data and Handle NAs
# Filter pheno to remove samples with missing Age or Sex
pheno <- pheno %>% 
  filter(!is.na(Age) & !is.na(Sex))

# Align edata columns to the filtered pheno rows
common_samples <- intersect(pheno$Sample, colnames(edata))
pheno <- pheno[match(common_samples, pheno$Sample), ]
edata <- edata[, match(common_samples, colnames(edata))]

# 2. Dynamic Model Construction
# Check if covariates have variance
has_age_var <- var(pheno$Age) > 0
# Check if Sex has more than 1 level present
has_sex_var <- length(unique(pheno$Sex)) > 1

# Build formula string based on valid covariates
covars <- c()
if (has_age_var) covars <- c(covars, "Age")
if (has_sex_var) covars <- c(covars, "Sex")

if (length(covars) > 0) {
  f_string <- paste("~", paste(covars, collapse = " + "))
  mod <- model.matrix(as.formula(f_string), data = pheno)
} else {
  # If no covariates vary, model is just intercept
  mod <- model.matrix(~ 1, data = pheno)
}

mod0 <- model.matrix(~ 1, data = pheno)

# 3. Check for Singularity (Collinearity)
# If Age and Sex are perfectly correlated, rank will be deficient
if (qr(mod)$rank < ncol(mod)) {
  warning("Design matrix is rank deficient. Reverting to intercept-only model.")
  mod <- model.matrix(~ 1, data = pheno)
}

# 4. Run SVA
# Only run SVA if we have covariates to protect (mod differs from mod0)
if (ncol(mod) > ncol(mod0)) {
  svobj <- sva(edata, mod, mod0, n.sv = 5)
  sv_factors <- as.data.frame(svobj$sv)
} else {
  # If no valid covariates exist, we cannot 'protect' them. 
  # Fill with zeros to prevent downstream errors.
  warning("No valid covariates to protect. Filling SVs with zeros.")
  sv_factors <- data.frame(matrix(0, nrow = nrow(pheno), ncol = 5))
}

colnames(sv_factors) <- paste0("SV", seq_len(ncol(sv_factors)))
pheno_with_svs <- cbind(pheno, sv_factors)

# Match samples with genotype data
matching_columns <- wgs_subset$V2[wgs_subset$V2 %in% colnames(edata)]
edata_subset <- edata[, matching_columns, drop = FALSE]

# PCA
pca_result <- prcomp(t(edata_subset), scale. = TRUE)
top_30_pcs <- pca_result$x[, seq_len(min(30, nrow(pca_result$x)))]
rownames(top_30_pcs) <- colnames(edata_subset)

# Create phenotype and covariate data
phenotype <- as.data.frame(t(edata_subset)) %>%
  rownames_to_column("Sample")

covs1 <- merge(
  metadata %>% filter(Sample %in% row.names(top_30_pcs)), 
  pheno_with_svs, 
  by = "Sample"
) %>%
  select(Sample, Sex.x, Age.x, SV1, SV2, SV3, SV4, SV5) %>%
  rename(Sex = Sex.x, Age = Age.x)

pcs <- as.data.frame(top_30_pcs) %>% rownames_to_column("Sample")
clusters <-  read.csv(cfile) %>%
  filter(X %in% colnames(edata)) %>%
  rename(Sample=X)

masterdf <- merge(merge(merge(merge(covs1, pcs, by = "Sample"), genotype_pcs, by = "Sample"), phenotype, by = "Sample") %>%
  mutate(Age_scaled = (Age - min(Age)) / (max(Age) - min(Age))), clusters, by="Sample")

valid_columns <- intersect(names(phenotype)[-1], names(masterdf))

phenotype <- masterdf %>%
  mutate(FID = 0, IID = Sample) %>%
  select(FID, IID, all_of(valid_columns)) %>%
  select(-FID) %>%
  column_to_rownames("IID") %>%
  t() %>%
  as.data.frame()

merged_data <- merge(
  phenotype %>% rownames_to_column("NAME"),
  gene_annotation,
  by = "NAME",
  all.x = TRUE
) %>%
  select(probe, chr, TSS, NAME, strand) %>%
  filter(!is.na(probe) & !is.na(chr) & !is.na(TSS))

# Update phenotype and match order
phenotype <- phenotype[rownames(phenotype) %in% merged_data$NAME, ]
phenotype <- phenotype[match(merged_data$NAME, rownames(phenotype)), ]

merged_data <- merged_data %>%
  filter(NAME %in% rownames(phenotype)) %>%
  arrange(match(NAME, rownames(phenotype)))

# Write outputs
write.table(as.data.frame(t(phenotype)) %>% rownames_to_column("Sample") %>%
              mutate(FID = 0, IID = Sample) %>% select(FID, IID, rownames(phenotype)),
            paste0(workdir, "/osca_input/Phenotype_", file, "_ocsa.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(merged_data %>% select(chr, NAME, TSS, probe, strand),
            paste0(workdir, "/osca_input/Upprobe_", file, ".opi"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(masterdf %>%
              mutate(FID = 0, IID = Sample) %>%
              select(FID, IID, Sex) %>%
              mutate(Sex = as.factor(Sex)),
            paste0(workdir, "/osca_input/cov1_", file, ".txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

cluster_names <- names(clusters)[-1]
covariates_ranked <- c("Age_scaled", "SV1", "G_PC1", "PC1", cluster_names,
                       "SV2", "G_PC2", "PC2", "SV3", "G_PC3", "PC3",
                       "SV4", "G_PC4", "PC4", "SV5", "G_PC5",
                       colnames(pcs)[6:ncol(pcs)])
# Reorder covariate columns into a priority list
masterdf <- masterdf %>%
  mutate(IID = Sample) %>%
  select(
    IID,
    any_of(covariates_ranked)
  )
# Remove constant columns so OSCA will work
const_cols <- sapply(masterdf, function(x) all(x == x[1]))
masterdf <- masterdf[!const_cols]
masterdf <- masterdf[seq_len(min(nrow(pcs) - 2, ncol(masterdf)))]

# Remove perfectly (or nearly) correlated columns so OSCA works
#cor_matrix <- cor(masterdf %>% select(!IID), use = "pairwise.complete.obs")
#high_cor <- which(abs(cor_matrix) > 0.999 & lower.tri(cor_matrix),
                  arr.ind = TRUE)
#masterdf <- masterdf %>% select(!rownames(high_cor))

write.table(add_column(masterdf, FID = 0, .before = 1),
            paste0(workdir, "/osca_input/cov2_", file, ".txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)
