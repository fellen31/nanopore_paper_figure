library(tidyverse)
library(DSS)
library(edgeR)
library(ggpubr)
library(gggenes)
library(MetBrewer)
library(cowplot)

# Read gene annotations
gff <- read_tsv(pipe("grep -v '#' /media/ssd_4tb/felix/projects/crebbp_methylation/promethion_methylation/gencode.v19.annotation.gff3"), col_names = c("chr", "V2", "V3", "start", "V5", "V6", "V7", "V8", "V9"))

# Filter for genes and get gene names 
gff <- gff %>% 
  filter(V3 == "gene") %>%
  dplyr::mutate(name = sub(".*gene_name=(.*?);.*", "\\1", V9))
  

# Read DSS-files   
file_paths <- list.files("/media/ssd_4tb/felix/projects/crebbp_methylation/promethion_methylation/autism_methylation/", pattern = "*.dss_16$", full.names = TRUE)

# Create a function to read tabular files without headers
read_tabular_file <- function(file_path) {
  read.table(file_path, col.names = c("chr", "pos", "N", "X"))
}

# Use lapply to read and store tabular files in a list
file_contents <- lapply(file_paths, read_tabular_file)

# Function to extract digits from filenames (use as sample names)
extract_digits <- function(string) {
  digits <- as.numeric(sub(".*op_\\d{3}_(\\d{3}).*", "\\1", string))
  formatted_digits <- sprintf("%03d", as.numeric(digits))
  formatted_digits
}

# Apply the function to each string in the list
sample_names <- lapply(file_paths, extract_digits)

# Display the extracted digits
print(sample_names)

# Make DDS object  
BSobj = makeBSseqData(file_contents, sampleNames = sample_names)

# Smoothing window is 500 bp by default...
make_dmr_plot <- function(arrangement, group_1, group_2, pval_t) {
  
  # Proband against all but father
  dmltest_004_against_all = DMLtest(BSobj, group1= c("004"), group2 = c("001", "003", "005", "006", "007", "008", "009", "010", "011", "012", "014", "015", "016", "017", "019", "020", "021"), 
                       smoothing=TRUE, ncores = 36)
  
  dmls_004_against_all = callDML(dmltest_004_against_all, p.threshold = pval_t)
  
  dmls_004_against_all %>% 
    as_tibble() %>% 
    dplyr::mutate(start = pos) %>%
    ggplot(aes(x = pos, y = -log10(pval), color = -log10(pval))) + 
    geom_point(alpha = .8) +
    theme_pubclean() -> dml_plot
  
  return(dml_plot)
}

make_gene_plot <- function (arrangement, gff_df) {
  gff_df %>%
    arrangement %>%
    tidytable::mutate(
      molecule = ifelse(V7 == "+", 1,2), 
      gene = name, 
      start = start, 
      end = V5, 
      strand = ifelse(V7 == "+", "forward", "reverse"), 
      orientation = ifelse(V7 == "-", 0, 1)
    ) %>%
    dplyr::select(gene, start, end, strand, orientation, molecule) %>% 
    as.data.frame() %>%
    ggplot(aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene, forward = orientation)) +
    geom_gene_arrow() +
    geom_gene_label(aes(label = gene)) +
    theme_genes() +
    theme(legend.position = "none") -> gene_plot
  
  return(gene_plot)
}

read_modkit_files <- function(path_to_files, name_pattern) {
  
  modkit_files = dir(path_to_files, pattern = name_pattern, recursive = TRUE)
  
  methylation <-tibble(file = modkit_files) %>%
    tidytable::mutate(file_contents = tidytable::map(file, 
                                                     ~ fread(paste0(path_to_files, .), sep = '\t', col.names = c("chr", "start", "end", "base", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18"))),
    ) %>%
    unnest(file_contents) %>%
    mutate(sample = sub(".*op_\\d{3}_(\\d{3}).*", "\\1", file))
  
  return(methylation)
}

methylation <- read_modkit_files("/media/ssd_4tb/felix/projects/crebbp_methylation/promethion_methylation/autism_methylation/", "*no_threshold.16$")

make_meth_plot <- function(methylation_df, arrangement, highlight_line, highlight_line2, n_cgs) {
  
  methylation_df %>%
    arrangement %>%
    tidytable::group_by(sample) %>%
    tidytable::arrange(start) %>%
    tidytable::mutate(meth = frollmean(X11, n_cgs, align = "center")) -> plot_data 
  
  plot_data %>% 
      ggplot(aes(start, meth, group = sample, color = sample)) +
      geom_vline(xintercept = 3688503, linetype = 2) + geom_vline(xintercept = 3899656, linetype = 2) + 
      geom_line(data=plot_data[which(sample==highlight_line)],linewidth = .8, alpha = 1) +
      geom_line(data=plot_data[which(sample==highlight_line2)],linewidth = .8, alpha = 1) +
      geom_line(alpha = .2, linewidth = 0.4) +
      ylab("Methylation (%)") +
      xlab("Position chr16") +
      theme_pubr() -> plot
  
  return(plot)
}

CREBBP <- function(df) {
  df %>%
    filter(chr == "chr16" | chr == "16") %>%
    filter(between(start, 3.*10^6, 4.2*10^6)) %>%
    arrange(start) -> df_filtered
  return(df_filtered)
}

p1 <- make_meth_plot(methylation, arrangement = CREBBP, "004", "002", 150) + theme_pubclean() + xlab("") + 
  scale_x_continuous(labels = scales::comma, limits = c(3688503-1*10^5,3899656+1*10^5)) + ylim(0,100) + theme(legend.position = "none")

p2 <- make_gene_plot(arrangement = CREBBP, gff %>% filter(name %in% c("DNASE1", "TRAP1", "CREBBP"))) + 
  scale_x_continuous(labels = scales::comma, limits = c(3688503-1*10^5,3899656+1*10^5)) + 
  theme_void() + 
  ylab("") + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none") + 
  ylim(1-1,2+1) 


p3 <- saved_dmr_plot + 
  scale_x_continuous(labels = scales::comma, limits = c(3688503-1*10^5,3899656+1*10^5)) + 
  ylim(0,45) +
  ylab(expression(-log[10](p))) +
  xlab("") + theme(legend.position = "none")

plot_grid(p1,p2,p3, ncol = 1, align = "v", axis = "left", rel_heights = c(6/13, 1/13, 6/13), labels = c('B', '', 'C'))

saved_dmr_plot <- make_dmr_plot(CREBBP, c("004"), c("002"), 0.01) 

sessionInfo()
