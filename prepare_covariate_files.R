args = commandArgs(trailingOnly=TRUE)
library(edgeR)
library(Rsubread)
library(Biobase)
library(preprocessCore)
library(ggplot2)
library(ggthemes)
library(glmnet)
library(reshape)
library(cowplot)
library(mvtnorm)
library(sva)
library(limma)

# Helper method to save DGE data structure to tab-deliminated text file
save_python_style_matrix <- function(counts, output_file, row_label_name) {
    #  Convert DGE data structure to matrix
    temp_mat <- as.matrix(counts)

    #  Edit colnames to include a header over the row-labels.
    revised_column_names <- colnames(temp_mat)
    revised_column_names[1] <- paste0(row_label_name,"\t",revised_column_names[1])

    write.table(temp_mat, output_file, quote=FALSE,col.names = revised_column_names, sep="\t")

}

# Write PC scores to output file
save_pcs <- function(sample_names, quant_expr, n, output_file) {
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first n pc's across all samples
    pc <- svd1$v[,1:n]
    colnames(pc) <- paste0("PC",1:n)
    rownames(pc) <- sample_names
    save_python_style_matrix(pc, output_file, "Sample_id")
}

save_sva <- function(sample_info, quant_expr, sva_output_file) {
    # Formatting changes required for SVA
    rownames(sample_info) <- sample_info$Sample_name
    colnames(quant_expr) <- sample_info$Sample_name

    # Model of interest (time)
    mod <-  model.matrix(~time,data=sample_info)
    mod0 <- model.matrix(~1,data=sample_info)

    # Estimate number of surrogate variables
    t.sv <- num.sv(quant_expr,mod)

    # Run sva
    svaobj <- sva(as.matrix(quant_expr),mod,mod0,n.sv=t.sv)

    # Extract factors. Matrix of dim (Num_samples)X(t.sv)
    sva_factors <- svaobj$sv

    # Add row and column labels
    rownames(sva_factors) <- sample_info$Sample_name
    colnames(sva_factors) <- paste0("factor",1:svaobj$n.sv)
    
    save_python_style_matrix(sva_factors, sva_output_file, "Sample_id")
}



add_covariate_column <- function(covariate_info, raw_covariates, column_num, column_name, column_type) {
    #  Number of samples
    N <- dim(covariate_info)[1]
    # If column is real valued
    if (column_type == "real_valued") {
        # Initialize column in covariate file
        covariate_info$my_var <- numeric(N)
        # Add name to column
        colnames(covariate_info)[length(colnames(covariate_info))] <- column_name

        temp_col_num <- length(colnames(covariate_info))
        # Loop through samples
        for (n in 1:N) {
            n_time <- covariate_info$time[n]
            n_cell_line <- covariate_info$cell_line[n]
            # Find corresponding row index in row_covariates file
            row_index <- which(raw_covariates$Line == n_cell_line & raw_covariates$Day == n_time)
            # Find value of this specific covariate
            covariate_info[n,temp_col_num] <- raw_covariates[row_index,column_num]
        }
    }
    if (column_type == "categorical") {
        # Determine discrete categories
        categories <- unique(factor(raw_covariates[,column_num]))
        num_categories <- length(categories)
        for (num_category in 1:num_categories) {
            category_name <- categories[num_category]
            new_column_name <- paste0(column_name,"_",num_category)

            # Initialize column in covariate file
            covariate_info$my_var <- numeric(N)
            # Add name to column
            colnames(covariate_info)[length(colnames(covariate_info))] <- new_column_name
            temp_col_num <- length(colnames(covariate_info))
            raw_covariates_temp <- (raw_covariates[,column_num] == category_name)*1.0
            for (n in 1:N) {
                n_time <- covariate_info$time[n]
                n_cell_line <- covariate_info$cell_line[n]
                # Find corresponding row index in row_covariates file
                row_index <- which(raw_covariates$Line == n_cell_line & raw_covariates$Day == n_time)
                # Find value of this specific covariate
                covariate_info[n,temp_col_num] <- raw_covariates_temp[row_index]
            }

        }

    }
    return(covariate_info)
}

add_covariate_column_categorical <- function(covariate_info, raw_covariates, column_num, column_name, column_type) {
    #  Number of samples
    N <- dim(covariate_info)[1]
    # If column is real valued
    # Initialize column in covariate file
    if (column_type == "categorical") {
        raw_covariates[,column_num] <- factor(raw_covariates[,column_num])
        covariate_info$my_var <- rep("null",N)
        covariate_info$my_var <- factor(covariate_info$my_var)
        levels(covariate_info$my_var) <- levels(raw_covariates[,column_num])
    }
    if (column_type == "real_valued") {
        covariate_info$my_var <- numeric(N)
    }
    # Add name to column
    colnames(covariate_info)[length(colnames(covariate_info))] <- column_name

    temp_col_num <- length(colnames(covariate_info))
    # Loop through samples
    for (n in 1:N) {
        n_time <- covariate_info$time[n]
        n_cell_line <- covariate_info$cell_line[n]
        if (n_cell_line == "18499") { # This was added to deal with the sample swap issue!! (becuase the covariate input data still contains the old sample)
            n_cell_line <- "19238"  
        }
        # Find corresponding row index in row_covariates file
        row_index <- which(raw_covariates$Line == n_cell_line & raw_covariates$Day == n_time)
        # Find value of this specific covariate
        covariate_info[n,temp_col_num] <- raw_covariates[row_index,column_num]
    }
    return(covariate_info)
}

add_multiqc_covariate <- function(covariate_info, multiqc, column_num, column_name) {
    N <- dim(covariate_info)[1]
    covariate_info$my_var <- numeric(N)

    # Add name to column
    colnames(covariate_info)[length(colnames(covariate_info))] <- column_name

    temp_names <- as.vector(covariate_info$Sample_name)
    temp_col_num <- length(colnames(covariate_info))
    for (n in 1:N) {
        sample_name <- temp_names[n]
        row_num <- which(multiqc$Sample == sample_name)
        covariate_info[n, temp_col_num] <- multiqc[row_num,column_num]
    }
    return(covariate_info)
}

# Binarize categorical variables
add_covariate_column2 <- function(covariate_info, column_num, column_name) {
    #  Number of samples
    N <- dim(covariate_info)[1]
    #  Vector of unique categores
    categories <- unique(factor(covariate_info[,column_num]))
    # Number of categories
    num_categories <- length(categories)
    # Create new column for each category (binary variable whether in category or not)
    for (num_category in 1:num_categories) {
        # Name of current category
        category_name <- categories[num_category]
        # What we want to name this column num
        new_column_name <- paste0(column_name,"_",num_category)

        # Initialize column in covariate file
        covariate_info$my_var <- numeric(N)
        # Add name to column
        colnames(covariate_info)[length(colnames(covariate_info))] <- new_column_name

        # Get the new column number
        temp_col_num <- length(colnames(covariate_info))
        # Fill in elements with 1 if sample belongs to category, 0 if not
        covariate_info[,temp_col_num] <- 1.0*(covariate_info[,column_num] == category_name)
    }
    return(covariate_info)
}


process_and_save_covariates_binary_categorical <- function(sample_info, raw_covariates, cov_output_file) {
    # Initialize covariate matrix to be a copy of sample_info
    covariate_info <- sample_info
    # Number of samples
    N <- dim(sample_info)[1]

    # Add columns (this part is fairly manual)
    covariate_info <- add_covariate_column(covariate_info, raw_covariates, 3, "RNA_concentration", "real_valued")
    covariate_info <- add_covariate_column(covariate_info, raw_covariates, 4, "RIN", "real_valued")
    covariate_info <- add_covariate_column(covariate_info, raw_covariates, 5, "RNA_extraction_day", "categorical")
    covariate_info <- add_covariate_column(covariate_info, raw_covariates, 6, "RNA_extraction_person", "categorical")
    covariate_info <- add_covariate_column(covariate_info, raw_covariates, 7, "RNA_extraction_round", "categorical")
    covariate_info <- add_covariate_column(covariate_info, raw_covariates, 8, "X1", "real_valued")
    covariate_info <- add_covariate_column(covariate_info, raw_covariates, 9, "H20", "real_valued")
    covariate_info <- add_covariate_column(covariate_info, raw_covariates, 10, "Differ_batch", "categorical")
    covariate_info <- add_covariate_column(covariate_info, raw_covariates, 11, "Lib_batch", "categorical")
    covariate_info <- add_covariate_column2(covariate_info, 8, "cellLine")


    # Remove un-informative columns
    covariate_info <-  cbind(covariate_info[,1],covariate_info[,7], covariate_info[,9:length(colnames(covariate_info))])
    # Add column name.. Got removed accidentally
    colnames(covariate_info)[1] <- "sample_names"
    colnames(covariate_info)[2] <- "time"


    save_python_style_matrix(covariate_info, cov_output_file, "Order_sample_num")
}

convert_from_date_to_days <- function(array) {
    N <- length(array)
    day <- numeric(length(array))
    for (n in 1:N) {
        element <- array[n]
        month <- as.numeric(strsplit(as.character(element),"/")[[1]][1])
        dayz <- as.numeric(strsplit(as.character(element),"/")[[1]][2])
        year <- strsplit(as.character(element),"/")[[1]][3]
        new_day <- 0
        if (is.na(year)) {
            year <- "17"
        }
        if (year == "17") {
            new_day <- new_day + 365
        }
        if (year == "2017") {
            new_day <- new_day + 366
        }
 
        new_day <- month*30 + dayz
        day[n] <- new_day
    }
    day <- day
    return(day)
}


process_and_save_covariates_categorical <- function(sample_info, raw_covariates, multiqc, cov_output_file) {
    # Initialize covariate matrix to be a copy of sample_info
    covariate_info <- sample_info
    # Number of samples
    N <- dim(sample_info)[1]


    # Format column labels
    multiqc$Sample <- as.vector(multiqc$Sample)
    multiqc$Sample <- substr(multiqc$Sample,1,nchar(multiqc$Sample)-7)


    # Add columns (this part is fairly manual)
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 3, "feeder_passage","real_valued")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 4, "feeder_free_passage","real_valued")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 5, "rna_extraction_conc","real_valued")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 6, "RIN","real_valued")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 7, "RNA_extraction_day","categorical")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 9, "RNA_extraction_person","categorical")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 11, "RNA_extraction_round","categorical")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 12, "volume_needed","real_valued")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 13, "water_added","real_valued")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 14, "differentiation_batch","categorical")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 16, "library_batch","categorical")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 22, "beating","real_valued")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 23, "line_beating","real_valued")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 28, "flash_freezing","real_valued")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 27, "cell_network_description","categorical")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 29, "sheets_vs_clumps_description","categorical")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 30, "debris_description","categorical")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 31, "cloudy_description","categorical")
    covariate_info <- add_covariate_column_categorical(covariate_info, raw_covariates, 32, "swirl_description","categorical")

    # Now add stats from multiqc
    covariate_info <- add_multiqc_covariate(covariate_info, multiqc, 2, "fastq_percent_duplicates")
    covariate_info <- add_multiqc_covariate(covariate_info, multiqc, 3, "percent_gc")
    covariate_info <- add_multiqc_covariate(covariate_info, multiqc, 4, "average_sequence_length")
    covariate_info <- add_multiqc_covariate(covariate_info, multiqc, 5, "total_sequences")
    covariate_info <- add_multiqc_covariate(covariate_info, multiqc, 6, "percent_fails")


    # Remove un-informative columns
    covariate_info <-  cbind(covariate_info[,1],covariate_info[,6], covariate_info[,5], covariate_info[,3], covariate_info[,7:length(colnames(covariate_info))])
    # Add column name.. Got removed accidentally
    colnames(covariate_info)[1] <- "sample_names"
    colnames(covariate_info)[2] <- "time"
    colnames(covariate_info)[3] <- "cell_line"
    colnames(covariate_info)[4] <- "library_size"


    save_python_style_matrix(covariate_info, cov_output_file, "Order_sample_num")
}




preprocess_total_expression_dir = args[1]  # Where total expression processed data 
metadata_input_file = args[2]  # Input file created by Katie/Reem that has all organized covariate information
covariate_dir = args[3]  # Output did to save covariate information
fastqc_dir = args[4]  # Input directory that contains fastqc/multiqc output files

#  Get sample information 
sample_info_file <- paste0(preprocess_total_expression_dir, "sample_info.txt")
sample_info <- read.table(sample_info_file, header=TRUE)

#  Get quantile normalized expression data
quantile_normalized_exp_file <- paste0(preprocess_total_expression_dir, "quantile_normalized.txt")
quant_expr <- read.csv(quantile_normalized_exp_file, header=TRUE, sep=" ")

#  Get rpkm expression_data
rpkm_exp_file <- paste0(preprocess_total_expression_dir, "rpkm.txt")
rpkm_expr <- read.csv(rpkm_exp_file, header=TRUE, sep=" ")

# Get multiqc general stats output file
multiqc_file <- paste0(fastqc_dir, "multiqc_data/multiqc_general_stats.txt")
multiqc <- read.table(multiqc_file,header=TRUE)

#  Get raw covariates
raw_covariates <- read.csv(metadata_input_file)


##############################################################################################################################
# Write PCs to output file
##############################################################################################################################
#  Number of pcs to save
n <- 10
#  Ouptut file to save PC loadings to 
pc_output_file <- paste0(covariate_dir, "principal_components_", n, ".txt")
save_pcs(sample_info$Sample_name, quant_expr, n, pc_output_file)


##############################################################################################################################
# Write surrogate variables via SVA to output file
##############################################################################################################################
#  Ouptut file to save SVA loadings to 
sva_output_file <- paste0(covariate_dir, "sva_loadings.txt")
save_sva(sample_info, quant_expr, sva_output_file)



##############################################################################################################################
# Write raw covariates  to output file
##############################################################################################################################
# cov_output_file <- paste0(covariate_dir, "processed_covariates_binary_categorical.txt")
# process_and_save_covariates_binary_categorical(sample_info, raw_covariates, cov_output_file)


cov_output_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
process_and_save_covariates_categorical(sample_info, raw_covariates, multiqc, cov_output_file)
