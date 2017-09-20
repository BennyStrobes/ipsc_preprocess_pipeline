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


# Visualize the percent of het snps that show biallelic expression. Color points by cell line
percent_biallelic_het_snps_scatter_cell_line <- function(ref_counts, total_counts, sample_info, percent_biallelic_ouptut_file, num_read_threshold) {
    N <- dim(sample_info)[1]  # Num samples
    #  Add new column to sample_info to keep track of the percent bi-allelic in each sample
    sample_info$percent_biallelic <- numeric(N)

    # Loop through each sample (compute percent biallleic in each)
    for (n in 1:N) {
        n_ref <- ref_counts[,n]  # reference allele counts for the nth sample
        n_total <- total_counts[,n]  # total counts for the nth sample
        # Consider only heterozygous sites (the nan check) and sites that have greater than num_read_threshold reads mapped totally
        n_observed_indices <- !is.nan(n_total) & (n_total > num_read_threshold)

        # Compute whether each sites is biallelic in terms of reads mapped
        biallelic_sites <- (n_ref[n_observed_indices] != n_total[n_observed_indices]) & (n_ref[n_observed_indices] != 0)

        # Compute percent of het snps that show biallelic expression
        sample_info$percent_biallelic[n] <- sum(biallelic_sites)/length(biallelic_sites)
    }

    # Now reorder percent_biallelic so they are ordered by cell line, and also time
    ordered_biallelic <- c()
    ordered_cell_lines <- c()
    unique_cell_lines <- unique(sample_info$cell_line)
    # Loop through each cell line
    for (i in 1:length(unique_cell_lines)) {
        i_cell_line <- unique_cell_lines[i]  # current cell line
        i_sample_info <- sample_info[sample_info$cell_line == i_cell_line,]  # Subset sample_info to only include samples belonging to this cell line

        i_times <- sort(i_sample_info$time)  # Now sort times within this cell line in numerical order
        #  Loop through sorted times
        for (j in 1:length(i_times)) {
            ij_time <- i_times[j]  # current time
            # calculate percent biallelic in this cell line at this time step
            ordered_biallelic <- c(ordered_biallelic, i_sample_info[i_sample_info$time == ij_time,]$percent_biallelic) 
            ordered_cell_lines <- c(ordered_cell_lines, i_cell_line)
        }

    }


    # Put data into data.frame for plotting
    df <- data.frame(sample_num = as.numeric(rownames(sample_info)), percent_biallelic = ordered_biallelic, cell_line = factor(ordered_cell_lines))

    #PLOT!
    scatter <- ggplot(df, aes(x = sample_num, y = percent_biallelic, colour = cell_line)) + geom_point() 
    scatter <- scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter  + labs(x = "Sample", y = "Percent Biallelic", colour = "Cell Line")
    ggsave(scatter, file=percent_biallelic_ouptut_file,width = 15,height=10.5,units="cm")
}



# Visualize the percent of het snps that show biallelic expression. Color points by cell line
percent_biallelic_het_snps_scatter_library_size <- function(ref_counts, total_counts, sample_info, percent_biallelic_ouptut_file, num_read_threshold) {
    N <- dim(sample_info)[1]  # Num samples
    #  Add new column to sample_info to keep track of the percent bi-allelic in each sample
    sample_info$percent_biallelic <- numeric(N)

    # Loop through each sample (compute percent biallleic in each)
    for (n in 1:N) {
        n_ref <- ref_counts[,n]  # reference allele counts for the nth sample
        n_total <- total_counts[,n]  # total counts for the nth sample
        # Consider only heterozygous sites (the nan check) and sites that have greater than num_read_threshold reads mapped totally
        n_observed_indices <- !is.nan(n_total) & (n_total > num_read_threshold)

        # Compute whether each sites is biallelic in terms of reads mapped
        biallelic_sites <- (n_ref[n_observed_indices] != n_total[n_observed_indices]) & (n_ref[n_observed_indices] != 0)

        # Compute percent of het snps that show biallelic expression
        sample_info$percent_biallelic[n] <- sum(biallelic_sites)/length(biallelic_sites)
    }

    # Now reorder percent_biallelic so they are ordered by cell line, and also time
    ordered_biallelic <- c()
    ordered_lib_sizes <- c()
    unique_cell_lines <- unique(sample_info$cell_line)
    # Loop through each cell line
    for (i in 1:length(unique_cell_lines)) {
        i_cell_line <- unique_cell_lines[i]  # current cell line
        i_sample_info <- sample_info[sample_info$cell_line == i_cell_line,]  # Subset sample_info to only include samples belonging to this cell line

        i_times <- sort(i_sample_info$time)  # Now sort times within this cell line in numerical order
        #  Loop through sorted times
        for (j in 1:length(i_times)) {
            ij_time <- i_times[j]  # current time
            # calculate percent biallelic in this cell line at this time step
            ordered_biallelic <- c(ordered_biallelic, i_sample_info[i_sample_info$time == ij_time,]$percent_biallelic) 
            ordered_lib_sizes <- c(ordered_lib_sizes, i_sample_info[i_sample_info$time == ij_time,]$lib.size)
        }

    }


    # Put data into data.frame for plotting
    df <- data.frame(sample_num = as.numeric(rownames(sample_info)), percent_biallelic = ordered_biallelic, library_size = ordered_lib_sizes)
    #PLOT!
    scatter <-  ggplot(df,aes(sample_num,percent_biallelic)) + geom_point(aes(colour=library_size)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    scatter <- scatter + scale_color_gradient(low="pink",high="blue") 
    scatter <- scatter + labs(x = "Sample", y = "Percent Biallelic", colour = "Lib. Size")
    ggsave(scatter, file=percent_biallelic_ouptut_file,width = 15,height=10.5,units="cm")

}

# Compute fraction of imputated genotype sites for each sample
fraction_of_hard_coded_genotype_sites <- function(het_prob_file, sample_info, output_file) {
    # Stream het prob file (each line is a site)
    stop = FALSE
    f = file(het_prob_file, "r")
    total_sites <- 0
    while(!stop) {
        next_line = readLines(f, n = 1)
        if(length(next_line) == 0) {  # Stopping criteria
            stop = TRUE
            close(f)
            break
        }
        data <- unlist(strsplit(next_line,"\t"))
        if (startsWith(data[1],'#CHROM')) {  # HEADER
            cell_lines <- data[10:length(data)]  # Keep track of names of cell lines
            hard_coded <- numeric(length(cell_lines))  # Keep track of how many of the genotype calls are hard-coded
        }
        if (startsWith(data[1],'#') == FALSE) {  # Standard line
            total_sites <- total_sites + 1 # Keep track of total sites
            probs <- as.numeric(data[10:length(data)])
            hard_coded_1 <- 1.0*(probs==1.0)  # Binary variable whether each cell line is hard coded 1 for this genotype
            hard_coded_0 <- 1.0*(probs==0.0)  # Binary variable whether each cell line is hard coded 0 for this genotype
            hard_coded <- hard_coded + hard_coded_1 + hard_coded_0  # Keep track
        }
    }

    df <- data.frame(cell_line = factor(cell_lines), fraction_hard_coded = hard_coded/total_sites)

    bar_plot <- ggplot(df, aes(x=cell_line, y=fraction_hard_coded, fill=cell_line)) + geom_bar(stat="identity")

    bar_plot <- bar_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90, vjust=.43,hjust = -1)) 
    bar_plot <- bar_plot + labs(colour="Cell Line",x = "Cell Line", y = "% hard-coded sites") + theme(legend.position="none")

    ggsave(bar_plot, file=output_file,width = 15,height=10.5,units="cm")
}

# Make boxplot of number of expressed het-snps per individual with one box for every read cutoff (that defines what is an expressed het-snp)
number_of_expressed_het_snps_per_individual <- function(total_counts, output_filer, het_thresh) {
    # Keep track of variables (initialization)
    num_expressed_sites <- c()
    thresholds <- c()
    cell_line <- c()
    
    #Get order of cel_line ids
    ordered_cell_lines <- substr(colnames(total_counts),2,6)

    # Compute number of expressed het_snps for each of these thresholds
    read_threshs <- c(2,6,10,15,30)
    for (iter in 1:length(read_threshs)) {
        read_thresh <- read_threshs[iter]
        # Compute which (site,samples) have expression
        binary_matrix <- total_counts > read_thresh
        # Create array of length number of samples that is the number of expressed sites in each sample
        num_expressed <- colSums(binary_matrix,na.rm=TRUE)

        # keep track of data
        num_expressed_sites <- c(num_expressed_sites, num_expressed)
        cell_line <- c(cell_line, ordered_cell_lines)
        thresholds <- c(thresholds, numeric(length(num_expressed)) + read_thresh)
    }
    # PLOT
    df <- data.frame(cell_line = factor(cell_line), num_expressed_sites = num_expressed_sites, thresholds=factor(thresholds))
    box_plot <- ggplot(df, aes(x=thresholds, y=num_expressed_sites)) + geom_boxplot(notch=TRUE)
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(x = "Read Threshold", y = "Expressed het sites / Sample", title=paste0("Het threshold = .",het_thresh))
    ggsave(box_plot, file=output_filer,width = 15,height=10.5,units="cm")
}

# Make boxplot of number of expressed het-snps per individual with one box for every read cutoff (that defines what is an expressed het-snp) and also per cell_line
number_of_expressed_het_snps_per_individual_cell_line_binned <- function(total_counts, output_filer, het_thresh) {
    # Keep track of variables (initialization)
    num_expressed_sites <- c()
    thresholds <- c()
    cell_line <- c()
    
    #Get order of cel_line ids
    ordered_cell_lines <- substr(colnames(total_counts),2,6)

    # Compute number of expressed het_snps for each of these thresholds
    read_threshs <- c(2,6,10)
    for (iter in 1:length(read_threshs)) {
        read_thresh <- read_threshs[iter]
        # Compute which (site,samples) have expression
        binary_matrix <- total_counts > read_thresh
        # Create array of length number of samples that is the number of expressed sites in each sample
        num_expressed <- colSums(binary_matrix,na.rm=TRUE)

        # keep track of data
        num_expressed_sites <- c(num_expressed_sites, num_expressed)
        cell_line <- c(cell_line, ordered_cell_lines)
        thresholds <- c(thresholds, numeric(length(num_expressed)) + read_thresh)
    }
    # PLOT
    df <- data.frame(cell_line = factor(cell_line), num_expressed_sites = num_expressed_sites, thresholds=factor(thresholds))
    box_plot <- ggplot(df, aes(x=thresholds, y=num_expressed_sites, fill=cell_line)) + geom_boxplot()
    box_plot <- box_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    box_plot <- box_plot + labs(fill= "Cell Line",x = "Read Threshold", y = "Expressed het sites / Sample", title=paste0("Het threshold = .",het_thresh))
    ggsave(box_plot, file=output_filer,width = 15,height=10.5,units="cm")
}


# Command line arguments
processed_allelic_counts_dir = args[1]  # Contains allelic counts
visualize_allelic_counts_dir = args[2]  # Directory to save output files to
genotype_dir = args[3]  # genotype file
preprocess_total_expression_dir = args[4]  # contains sample_info file


#  Get sample information 
sample_info_file <- paste0(preprocess_total_expression_dir, "sample_info.txt")
sample_info <- read.table(sample_info_file, header=TRUE)


het_thresh <- 999 #  Heterozygous probability threshold

# Load in data 
# total_counts_file contains total number of reads mapping to each site in each sample
# ref_counts_file contains number of reads mapping to reference allele at each site in each sample
total_counts_file <- paste0(processed_allelic_counts_dir, "allelic_counts_gene_mapped_het_prob_", het_thresh, "_total_counts.txt")
ref_counts_file <- paste0(processed_allelic_counts_dir, "allelic_counts_gene_mapped_het_prob_", het_thresh, "_ref_counts.txt")
total_counts <- read.table(total_counts_file, header=TRUE, sep="\t")
total_counts <- total_counts[,2:dim(total_counts)[2]]  # Skip first column due to transition from python (no information lost)
ref_counts <- read.table(ref_counts_file, header=TRUE, sep="\t")
ref_counts <- ref_counts[,2:dim(ref_counts)[2]] # Skip first column due to transition from python (no information lost)


# Compute fraction of imputated genotype sites for each sample
het_prob_file <- paste0(genotype_dir, "YRI_het_prob_genotype.vcf")
fraction_hard_coded_output_file <- paste0(visualize_allelic_counts_dir, "percent_hard_coded_genotypes_by_cell_line.png")
#fraction_of_hard_coded_genotype_sites(het_prob_file, sample_info, fraction_hard_coded_output_file)


# Visualize the percent of het snps that show biallelic expression. Color points by cell line
num_read_threshold <- 5  # Only consider sites that have at least num_read_threshold reads mapping to both alleles
percent_biallelic_ouptut_file <- paste0(visualize_allelic_counts_dir, "percent_biallelic_het_snps_by_cell_line_",het_thresh, ".png")
#percent_biallelic_het_snps_scatter_cell_line(ref_counts, total_counts, sample_info, percent_biallelic_ouptut_file, num_read_threshold)

# Visualize the percent of het snps that show biallelic expression. Color points by total read-depth
num_read_threshold <- 5  # Only consider sites that have at least num_read_threshold reads mapping to both alleles
percent_biallelic_ouptut_file <- paste0(visualize_allelic_counts_dir, "percent_biallelic_het_snps_by_library_size_",het_thresh, ".png")
#percent_biallelic_het_snps_scatter_library_size(ref_counts, total_counts, sample_info, percent_biallelic_ouptut_file, num_read_threshold)




# Make boxplot of number of expressed het-snps per individual with one box for every read cutoff (that defines what is an expressed het-snp)
number_of_expressed_het_snps_output_file <- paste0(visualize_allelic_counts_dir, "number_of_expressed_het_snps_per_individual_boxplot_",het_thresh,".png")
#number_of_expressed_het_snps_per_individual(total_counts, number_of_expressed_het_snps_output_file, het_thresh)

# Make boxplot of number of expressed het-snps per individual with one box for every read cutoff (that defines what is an expressed het-snp) and also per cell_line
number_of_expressed_het_snps_output_file <- paste0(visualize_allelic_counts_dir, "number_of_expressed_het_snps_per_individual_cell_line_boxplot_",het_thresh,".png")
number_of_expressed_het_snps_per_individual_cell_line_binned(total_counts, number_of_expressed_het_snps_output_file, het_thresh)

