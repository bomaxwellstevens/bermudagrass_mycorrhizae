library(microeco)
library(devtools)
library(qiime2R)
library(magrittr)
library(ggplot2)
library(picante)
library(agricolae)
theme_set(theme_bw())

qiimed2meco <- function(ASV_data, sample_data, taxonomy_data, phylo_tree = NULL){
  # Read ASV data
  ASV <- as.data.frame(read_qza(ASV_data)$data)
  #  Read metadata
  metadata <- read_q2metadata(sample_data)
  rownames(metadata) <- as.character(metadata[, 1])
  # Read taxonomy table
  taxa_table <- read_qza(taxonomy_data)
  taxa_table <- parse_taxonomy(taxa_table$data)
  # Make the taxonomic table clean, this is very important.
  taxa_table %<>% tidy_taxonomy
  # Read phylo tree
  if(!is.null(phylo_tree)){
    phylo_tree <- read_qza(phylo_tree)$data
  }
  dataset <- microtable$new(sample_table = metadata, tax_table = taxa_table, otu_table = ASV, phylo_tree = phylo_tree)
  dataset
}
setwd('C:/Users/bo.stevens/OneDrive - USDA/2022/OSU Collab/QIIME2.2')

meco_dataset_16S <- qiimed2meco(ASV_data = "table-paired.qza", sample_data = "mapping_file/mapping_file_types.txt", taxonomy_data = "taxonomy-paired.qza", phylo_tree = "rooted-tree-paired.qza")
meco_dataset_16S$tidy_dataset()
meco_dataset_16S$tax_table[1:5, 1:3]
meco_dataset_16S$sample_sums() %>% range
# Could rarefy here
meco_dataset_16S$cal_abund()

dir.create("alpha_diversity")
meco_dataset_16S$cal_alphadiv(PD = TRUE)
meco_dataset_16S$save_alphadiv(dirpath = "alpha_diversity")

t16S <- trans_abund$new(dataset = meco_dataset_16S, taxrank = "Genus", ntaxa = 10)
t16S$plot_bar(others_color = "grey70", facet = "genotype_1", xtext_keep = FALSE, legend_text_italic = FALSE)

t16S <- trans_abund$new(dataset = meco_dataset_16S, taxrank = "Order", ntaxa = 12, groupmean = "genotype_1")
t16S$plot_bar(others_color = "grey70", legend_text_italic = FALSE)

meco_dataset_16S$cal_alphadiv(PD = TRUE)
t16S<- trans_alpha$new(dataset = meco_dataset_16S, group = "genotype_1")
t16S$alpha_stat
write.csv(t16S$alpha_stat, "alpha_diversity/alpha_stats.csv")
t16S$cal_diff(method = "anova")
# return t1$res_alpha_diff
t16S$res_alpha_diff
write.csv(t16S$res_alpha_diff, "alpha_diversity/alpha_diff_stats.csv")

dataset16S <- meco_dataset_16S$merge_samples(use_group = "genotype_1")
t16S <- trans_venn$new(dataset16S, ratio = "seqratio")
t16S$plot_venn(petal_plot = TRUE)
#LDA
t16S <- trans_diff$new(dataset = meco_dataset_16S, method = "lefse", group = "genotype_1", alpha = 0.01, lefse_subgroup = NULL)
t16S$plot_lefse_bar(LDA_score = 4)
#

