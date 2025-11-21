library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(qs)
library(DT)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
source("R/functions.R")
source("R/plotPCA2.R")
source("R/volcano.R")
data_dir <- "./data/processed_data/counts_data" ## the current folder
matrix_file <- "GSE234563_Matrix.txt.gz"

mat <- read.table(fs::path(data_dir,
                         matrix_file),
                header = TRUE,
                sep = "\t",
                quote = ""
)
dim(mat)

dt <- DT::datatable(mat |> head(),options(scrollX=TRUE))
dt

##First of all, we **convert the column names to all-lowercase**, to avoid case issues:
colnames(mat) <- tolower(colnames(mat))

colData <- strsplit(colnames(mat),"_")
head(colData, n=2)


colData <- Reduce(rbind,colData)
head(colData, n=2)

## Setting Up the Input for DESeq2
# Step 3: We know what the row and column names should be, so let's set them:
dimnames(colData) <- list(colnames(mat),
                          c("age","sex","genotype","tissue","replicate")
)
head(colData, n=2)

# Step 4: The DESeq2 colData input should be a data.frame:
colData <- as.data.frame(colData)


colData$group <-
  purrr::reduce(colData[,c("age","sex","genotype","tissue")],
                paste, sep=".")
apply(colData[,c("age","sex","genotype","tissue")],1,paste,collapse=".")

# Show a bit more of the final result:
head(colData, n=9)

tissue_of_interest <- "intestine"

keep <- colData$tissue == tissue_of_interest

# Subset the counts and metadata simultaneously
mat_intestine       <- mat[, keep]
colData_intestine <- colData[keep, , drop = FALSE]

# Drop unused factor levels
colData_intestine <- colData_intestine |>
  mutate(
    age      = factor(age),
    sex      = factor(sex),
    genotype = factor(genotype),
    group    = factor(group),
    tissue   = droplevels(factor(tissue)),
    replicate = factor(replicate),
    age_sex = interaction(sex, age)
  )

table(colData_intestine$tissue)

##Then we **discard genes with a total count < 10**, as was done in the paper:
mat_intestine <- ceiling(mat_intestine)
min_count <- 10
mat_intestine <- mat_intestine[rowSums(mat_intestine) >= min_count,]

nrow(mat_intestine)
##DeSeq2 Analysis: Modeling the Data With Genotype Only
dds_int <-
  # High level step 1: set up the object structure
  DESeqDataSetFromMatrix(countData = mat_intestine,
                         colData = colData_intestine,
                         design= ~ 1 + genotype + sex + age
  ) |>
  # High level step 2: "do the math" (--> Ali!)
  #                    size factors factors,
  #                    dispersion estimates,
  #                    statistical tests ...
  DESeq()

##At this point we can have the available test results listed:
resultsNames(dds_int)

## A Bold Try: Model the Data With Genotype Only
##We want to see the genotype effect ...

res_int_genotype <-
  results(dds_int,
          name="genotype_nmrhas2_vs_creer"
  )

res_int_sex <-
  results(dds_int,
          name="sex_male_vs_female"
  )

res_int_age <-
  results(dds_int,
          name="age_young_vs_old"
  )

saveRDS(res_int_genotype, file="res_int_genotype.rds")
saveRDS(res_int_sex, file="res_int_sex.rds")
saveRDS(res_int_age, file="res_int_age.rds")

dt <- DT::datatable(
  res_int_genotype |> as.data.frame() |>
    dplyr::filter(!is.na(padj), padj <= 0.05) |>
    dplyr::arrange(padj),

  options=list(scrollY="500px")
)
dt <- dt |>
  DT::formatSignif(
    purrr::map_lgl(dt$x$data, is.numeric),
    digits=2
  )

## A Bold Try: Model the Data With Genotype Onl
dt


dt_age <- DT::datatable(
  res_int_age |> as.data.frame() |>
    dplyr::filter(!is.na(padj), padj <= 0.05) |>
    dplyr::arrange(padj),

  options=list(scrollY="500px")
)

dt_age <- dt_age |>
  DT::formatSignif(
    purrr::map_lgl(dt_age$x$data, is.numeric),
    digits=2
  )

## A Bold Try: Model the Data With Genotype Onl
dt_age

# convert numeric columns to scientific notation

##boxplots of raw and normalized counts for all samples
library(dplyr)
library(reshape2)
library(ggplot2)

df <-
  rbind(reshape2::melt(log10(counts(dds_int,normalized=FALSE)+1)) |>
          mutate(mode="raw"),

        reshape2::melt(log10(counts(dds_int,normalized=TRUE)+1)) |>
          mutate(mode="normalized")
  ) |>
  # rename to make variable names more meaningful:
  dplyr::rename(gene=Var1, sample=Var2, count=value) |>

  # remove the trailing sample ID to create a common ID for all replicates
  # of a given factor level combination:
  mutate(group=sub("_\\d+$","",sample)) |>

  # ggplot needs to know the order of sub-plots ("raw" should come first)
  mutate(mode = factor(mode,
                       levels=c("raw","normalized")
  )
  )

## Then plot it:

ggplot(df,
       aes(y=sample,x=count,fill=group)
) +
  geom_boxplot(show.legend=FALSE) +
  facet_wrap(vars(mode),ncol=2) +
  theme_classic() +
  theme(axis.text=element_text(size=2))

##PCA and MA-PLots
plotMA(res_int_genotype, ylim = c(-5, 5))

##PCA
vsd <- vst(dds_int,)
saveRDS(vsd,file="vsd.rds")

## What Drives the Structure?
#The unlabelled dataset shows some very clear structures:
plotPCA2(vsd,
         intgroup = "group",
         color = group,
         shape = sex,
         ntop=500,
         PCs=c(x=1,y=2))


# Convert results to a data frame (important!)
res_df <- as.data.frame(res_int_genotype)

# Save to CSV
write.csv(
  res_df,
  file = "results_genotype_intestine_DESeq2.csv",
  row.names = TRUE
)

##volcano plots

library(dplyr)
library(tibble)

# Prepare DESeq2 results genotype
df_volcano <- res_int_genotype |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  dplyr::select(gene, log2FoldChange, padj)

# Plot volcano
volcano(
  df = df_volcano,
  use_cols = c(1, 2, 3),     # gene, log2 FC, p-value
  lfc_thrsh = 1,
  p_thrsh = 0.05,
  label.show = T,
  plotTheme = 0
)

# Prepare DESeq2 results sex
df_volcano <- res_int_sex |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  dplyr::select(gene, log2FoldChange, padj)

# Plot volcano
volcano(
  df = df_volcano,
  use_cols = c(1, 2, 3),     # gene, log2 FC, p-value
  lfc_thrsh = 1,
  p_thrsh = 0.05,
  label.show = T,
  plotTheme = 0
)

# Prepare DESeq2 results age
df_volcano <- res_int_age |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  dplyr::select(gene, log2FoldChange, padj)

# Plot volcano
volcano(
  df = df_volcano,
  use_cols = c(1, 2, 3),     # gene, log2 FC, p-value
  lfc_thrsh = 1,
  p_thrsh = 0.05,
  label.show = T,
  plotTheme = 0
)
##ORA and GSEA enrichment

uni_genes <- results(dds_int) |>
  as.data.frame() |>
  filter(!is.na(padj)) |>
  rownames()
de_genes_up <- results(dds_int) |>
  as.data.frame() |>
  filter(!is.na(padj), padj <= 0.05, log2FoldChange > 0) |>
  rownames()
de_genes_dn <- results(dds_int) |>
  as.data.frame() |>
  filter(!is.na(padj), padj <= 0.05, log2FoldChange < 0) |>
  rownames()

ego_result_up <- enrichGO(
  gene = de_genes_up,
  universe = uni_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

DT::datatable(head(ego_result_up))

ego_result_dn <- enrichGO(
  gene = de_genes_dn,
  universe = uni_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

DT::datatable(ego_result_dn|>data.frame())

# Visualizations
# x="Count"/"GeneRatio", color="p.adjust"/"qvalue", size="GeneRatio"/"Count"

dotplot(ego_result_up, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Upregulated Genes : Dotplot")
dotplot(ego_result_dn, showCategory=20) +
  ggtitle("GO Enrichment Analysis - Downregulated Genes : Dotplot")

barplot(ego_result_up, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Upregulated Genes : Barplot")
barplot(ego_result_dn, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Downregulated Genes : Barplot")

cnetplot(ego_result_up, categorySize="pvalue", foldChange=NULL) +
  ggtitle("GO Enrichment Analysis - Upregulated Genes : Cnetplot")
cnetplot(ego_result_dn, categorySize="pvalue", foldChange=NULL) +
  ggtitle("GO Enrichment Analysis - Downregulated Genes : Cnetplot")

ego_result_up <- enrichplot::pairwise_termsim(ego_result_up)
ego_result_dn <- enrichplot::pairwise_termsim(ego_result_dn)
emapplot(ego_result_up, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Upregulated Genes : Emapplot")
emapplot(ego_result_dn, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Downregulated Genes : Emapplot")


result_up <- results(dds_int) |>
  as.data.frame() |>
  filter(!is.na(padj), padj <= 0.05, log2FoldChange > 0) |>
  select(log2FoldChange) |>
  tibble::rownames_to_column(var = "gene") |>
  pull(log2FoldChange,gene)

cnetplot(ego_result_up, foldChange = result_up) +
  ggtitle("GO Enrichment Analysis - Upregulated Genes : Cnetplot with log2FC")

##GSEA

res_int.genotype <- results(dds_int, name="genotype_nmrhas2_vs_creer")


I <- (res_int.genotype |> as.data.frame() |>
        filter(!is.na(padj)) |> arrange(desc(stat)) |>
        as.matrix())[,"stat"]

gse_int_genotype_0.05 <- runGSEA(I,org.Mm.eg.db,"SYMBOL") ## BP NES+ fatty acid and organic acid catabolic processes

dotplot(gse_int_genotype_0.05$MF, showCategory = 20) +
  ggtitle("GSEA – top enriched BP pathways")

ridgeplot(gse_int_genotype_0.05$BP, showCategory = 20) +
  ggtitle("GSEA – ridgeplot of enriched pathways")


str(gse_int_genotype_0.05, max.level = 1)

colData(dds)$tissue |> table()

renv::snapshot()
