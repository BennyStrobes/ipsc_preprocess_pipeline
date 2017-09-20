args = commandArgs(trailingOnly=TRUE)
library(edgeR)
library(Rsubread)
library(Biobase)
library(preprocessCore)
library(ggplot2)
library(ggthemes)
library(reshape)
library(cowplot)
library(mvtnorm)




ref_file = args[1]
tot_file = args[2]
output_dir = args[3]
rna_seq_sample_id = args[4]

total_counts <- read.table(tot_file, header=TRUE, sep="\t")
total_counts <- total_counts[,2:dim(total_counts)[2]]  # Skip first column due to transition from python (no information lost)

ref_counts <- read.table(ref_file, header=TRUE, sep="\t")
ref_counts <- ref_counts[,2:dim(ref_counts)[2]]


N <- dim(total_counts)[2]

percent_biallelic <- numeric(N)
num_read_threshold <- 5
for (n in 1:N) {
    n_ref <- ref_counts[,n]
    n_total <- total_counts[,n]
    n_observed_indices <- !is.nan(n_total) & (n_total > num_read_threshold)

    # Compute whether each sites is biallelic in terms of reads mapped
    biallelic_sites <- (n_ref[n_observed_indices] != n_total[n_observed_indices]) & (n_ref[n_observed_indices] != 0)

    # Compute percent of het snps that show biallelic expression
    percent_biallelic[n] <- sum(biallelic_sites)/length(biallelic_sites)
}

cell_lines <- substring(colnames(total_counts)[1:length(colnames(total_counts))], 2)
hit <- percent_biallelic > .8
# Put data into data.frame for plotting
df <- data.frame(sample_num = 1:length(cell_lines), percent_biallelic = percent_biallelic)

#PLOT!
scatter <- ggplot(df, aes(x = sample_num, y = percent_biallelic)) + geom_point() 
scatter <- scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
scatter <- scatter  + labs(x = "Line", y = "Percent Biallelic", title=paste0("Sample=",rna_seq_sample_id," / matched Geno=",cell_lines[hit]))
ggsave(scatter, file=paste0(output_dir,rna_seq_sample_id,"_percent_biallelic.png"),width = 15,height=10.5,units="cm")
