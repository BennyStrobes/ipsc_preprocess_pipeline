args = commandArgs(trailingOnly=TRUE)
library(edgeR)
library(Rsubread)
library(Biobase)
library(preprocessCore)
library(ggplot2)
library(ggthemes)
library(glmnet)
library(reshape)
library(rstan)
library(cowplot)
library(mvtnorm)


BETABINOMIAL_GLM=stan_model(file="sparse_betabinomial_glm.stan", save_dso=T, auto_write=T)

#  Plot first two PC's. Color points by time step
plot_pca_time_step <- function(sample_info, quant_expr, output_file) {

    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,1]
    pc2 <- svd1$v[,2]

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, time_step = sample_info$time)

    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue")


    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")

}


plot_library_size <- function(sample_info, library_size_output_file) {
    #  Get unique cell lines
    cell_lines <- sort(unique(sample_info$cell_line))
    num_cell_lines <- length(cell_lines)

    #  Initialize arrays to store information
    ordered_library_sizes <- c()
    ordered_cell_lines <- c()
    ordered_time_steps <- c()
    ordered_sample_name <- c()
    
    # Loop through cell lines and samples
    for (cell_line_iter in 1:num_cell_lines) {
        for (time_step in 0:15) {
            # Extract cell line of ith iteration
            cell_line <- cell_lines[cell_line_iter]
            # Find which row this corresponds to
            correct_row <- which(sample_info$cell_line == cell_line & sample_info$time == time_step)
            # Check to make sure cellLine_timeStep exists in our data
            if (length(correct_row) > 0) {
                # If it does, append arrays
                ordered_library_sizes <- c(ordered_library_sizes,sample_info$lib.size[correct_row])
                ordered_cell_lines <- c(ordered_cell_lines, sample_info$cell_line[correct_row])
                ordered_time_steps <- c(ordered_time_steps, sample_info$time[correct_row])
                ordered_sample_name <- c(ordered_sample_name, sample_info$Sample_name[correct_row])
            }
        }
    }
    df <- data.frame(library_size = ordered_library_sizes, cell_line = factor(ordered_cell_lines), time_step = ordered_time_steps, sample_name = ordered_sample_name)

    bar_plot <- ggplot(df, aes(x=sample_name, y=library_size, fill=cell_line)) + geom_bar(stat="identity")

    bar_plot <- bar_plot + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    bar_plot <- bar_plot + labs(colour="Cell Line",x = "Sample", y = "Library Size")

    ggsave(bar_plot, file=library_size_output_file,width = 15,height=10.5,units="cm")


}


plot_pca_real_valued_gene_filled <- function(sample_info, quant_expr, ensamble_id, gene_name, pc_num1, pc_num2, output_file) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,pc_num1]
    pc2 <- svd1$v[,pc_num2]

    row_label <- which(rownames(quant_expr) == ensamble_id)
    quant_expr <- as.matrix(quant_expr)

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, time_step = as.vector(quant_expr[row_label,]))


    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue") + labs(colour="Expression",x = paste0("PC",pc_num1), title = gene_name,y = paste0("PC",pc_num2))


    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")
}


gene_time_course_line_plot_grouped_by_cell_line <- function(sample_info, quant_expr, ensamble_id, gene_name, line_plot_file) {
    row_label <- which(rownames(quant_expr) == ensamble_id)
    quant_expr <- as.matrix(quant_expr)

    df <- data.frame(time = sample_info$time, expression = as.vector(quant_expr[row_label,]), cell_line = factor(sample_info$cell_line))

    line_plot <- ggplot(df, aes(x=time, y=expression, group=cell_line)) + geom_line(aes(color=cell_line)) +
                geom_point(aes(color=cell_line)) +
                theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                labs(colour="Cell Line",x = "Time Step", y = "Normalized Expression", title = gene_name)
    ggsave(line_plot, file=line_plot_file,width = 15,height=10.5,units="cm")

}

pc_gene_scatter <- function(sample_info, quant_expr, ensamble_id, gene_name, time_step, pc_num, pc_gene_scatter_output_file) {
    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of pcs of interest
    pc_scores <- svd1$v[,pc_num]

    # Row corresponding to gene of interest
    row_label <- which(rownames(quant_expr) == ensamble_id)
    # Get matrix into correct format
    quant_expr <- as.matrix(quant_expr)


    df <- data.frame(time = sample_info$time, expression = as.vector(quant_expr[row_label,]), pc_scores = pc_scores, cell_line = factor(sample_info$cell_line))
    
    df <- df[df$time == time_step,]

    #PLOT!
    pca_scatter <- ggplot(df, aes(x = expression, y = pc_scores, colour = cell_line)) + geom_point() 
    pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(colour="Cell Line",x = paste0(gene_name, " expression"), y = paste0("PC",pc_num), title = paste0("Time step ", time_step))
    ggsave(pca_scatter, file=pc_gene_scatter_output_file,width = 15,height=10.5,units="cm")

}



#  Plot first two PC's. Color points by cell_line
plot_pca_categorical_covariate <- function(sample_info, quant_expr, output_file, covariate, covariate_name,pc_num1,pc_num2) {

    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,pc_num1]
    pc2 <- svd1$v[,pc_num2]

    # Put all information into data structure
    df <- data.frame(pc1 = pc1, pc2 = pc2, covariate = covariate)

    #PLOT!
    pca_scatter <- ggplot(df, aes(x = pc1, y = pc2, colour = covariate)) + geom_point() 
    pca_scatter <- pca_scatter + theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + labs(colour=covariate_name,x = paste0("PC",pc_num1), y = paste0("PC",pc_num2))
    ggsave(pca_scatter, file=output_file,width = 15,height=10.5,units="cm")

}

# Make plot showing variance explained of first n pcs
plot_pca_variance_explained <- function(sample_info, quant_expr, n, output_file) {
    sv <- svd(as.matrix(quant_expr))
    pdf(output_file)
    var_expl_plot <- plot(sv$d^2/sum(sv$d^2), xlim = c(0, n), type = "b", pch = n+1, xlab = "principal components", ylab = "variance explained")
    print(var_expl_plot)
    dev.off()
}





# Helper method to make pca plot
make_pca_plot <- function(i, cell_lines, pc1, pc2, sample_info) {
    # Get name of ith cell line
    i_cell_line <- cell_lines[i]

    # Get indices  of all samples that are from ith cell line
    i_indices <- sample_info$cell_line == i_cell_line

    #  Extract 1st two pcs of all samples that belong to the ith cell lein
    i_pc1 <- pc1[i_indices]
    i_pc2 <- pc2[i_indices]

    # Get time steps of these samples
    i_time <- sample_info$time[i_indices]

    # Put into compact data frame for plotting
    df <- data.frame(pc1 = i_pc1, pc2 = i_pc2, time_step = i_time)

    #PLOT!
    pca_scatter <-  ggplot(df,aes(pc1,pc2)) + geom_point(aes(colour=time_step)) + theme(text = element_text(size=12), axis.text = element_text(size=8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    pca_scatter <- pca_scatter + scale_color_gradient(low="pink",high="blue") + ggtitle(i_cell_line) + xlim(min(pc1)-.01,max(pc1)+.01) + ylim(min(pc2)-.01,max(pc2) + .01)

    return(pca_scatter)
}


#  Perform PCA. Make one plot for each cell line of first two pcs.
plot_pca_seperate_cell_lines <- function(sample_info, quant_expr, output_file,pc_num1,pc_num2) {

    #  Compute singular value decomposition
    svd1 <- svd(as.matrix(quant_expr))

    #  Scores of first 2 pc's across all samples
    pc1 <- svd1$v[,pc_num1]
    pc2 <- svd1$v[,pc_num2]


    #  Get unique cell lines
    cell_lines <- sort(unique(sample_info$cell_line))
    num_cell_lines <- length(cell_lines)


    
    #  Make pc plot for each cell line seperately (not automated yet..
     
    p1 <- make_pca_plot(1, cell_lines, pc1, pc2, sample_info) + labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) # cell line 1
    p2 <- make_pca_plot(2, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) + theme(legend.position="none") # cell line 2..
    p3 <- make_pca_plot(3, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p4 <- make_pca_plot(4, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p5 <- make_pca_plot(5, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p6 <- make_pca_plot(6, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p7 <- make_pca_plot(7, cell_lines, pc1, pc2, sample_info)+labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) + theme(legend.position="none")
    p8 <- make_pca_plot(8, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p9 <- make_pca_plot(9, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
    p10 <- make_pca_plot(10, cell_lines, pc1, pc2, sample_info)+ labs(x=paste0("PC",pc_num1),y=paste0("PC",pc_num2)) +theme(legend.position="none")
 

    legend <- get_legend(p1)

    # Merge all cell lines into one plot
    pdf(output_file)
    gg <- plot_grid(p1 + theme(legend.position="none"),p2,p3,p4,p5,p6,p7,p8,p9,p10,nrow=4,ncol=3)
    combined_gg <- ggdraw() + draw_plot(gg,0,0,1,1) + draw_plot(legend,.83,0,1,.3)
    print(combined_gg)
    dev.off()
}



# Fit the current model using training data
fit_model <- function(x_train, y_train, lambda, model_type) {
    if (model_type == "linear_regression") {
        #  Fit lasso linear regression with glmnet
        fit <- glmnet(x=x_train, y=y_train, alpha=1, lambda=lambda)
        #  Extract coefficient matrix of dim (num_genes + 1, 1)
        beta <- as.matrix(fit$beta)
        #  Include intercept..
        beta <- rbind(fit$a0,beta)
    } else if (model_type == "beta_binomial_regression") {
        #  Organize data for stan
        concShape <- 1.0001  # Diffuse prior hyperparameters for BB
        concRate <- 1e-4  # Diffuse prior hyperparameters for BB
        ns <- numeric(length(y_train)) + 15  # n in beta binomial is max of all time steps
        dat <- list(N=length(y_train), P=ncol(x_train), ys=y_train, x=x_train, concShape=concShape, concRate=concRate, lambda=lambda, ns=ns)

        ############################################
        # Propper initialization
        ###########################################
        rat <- dat$ys/dat$ns # ratios
        # moment estimator of concentration parameter
        conc <- pmin( 1/var(rat, na.rm=T), 1000 )
        m <- pmin( pmax( mean(rat, na.rm=T), 1/1000 ), 1-1/1000)
        betaInit <- numeric(ncol(x_train))
        betaInit[1] <- log(m/(1.0-m))
        init=list(conc=conc, beta=as.array(betaInit),alpha=0)


        #  Run stan second order optimization
        fit <- optimizing(BETABINOMIAL_GLM, data=dat,iter=10000, init=init, algorithm="LBFGS", hessian=T, as_vector=F, verbose=TRUE)

        #  Extract Parameters
        beta <- as.matrix(fit$par$beta)
        #  Include intercept..
        beta <- rbind(fit$par$alpha, beta)
    }
    return(beta)
}

#  Use fitted model to make predictions on testing data
model_prediction <- function(x_test, model_type, beta) {
    #  Add intercept to model
    column_of_ones <- (numeric(nrow(x_test)) +1)
    x_test<- cbind(column_of_ones, x_test)
    if (model_type == "linear_regression") {
        # Compute linear predictions
        y_predicted <- x_test %*% beta
    } else if (model_type == "beta_binomial_regression") {
        y_predicted <- (1/(1 + exp(-(x_test %*% beta))))*15
    }

    # Return in vector form
    return(y_predicted[,1])
}

# Determine optimal choice of lambda through k-fold cross validation
# x is design matrix
# y is response vector
# lambdas is 1d grid of possible lambdas
# k puts the k in k-fold
# model_type is type of glm
select_lambda_through_k_fold_cv <- function(x, y, lambdas, k, model_type) {
    #  Shuffle data
    shuffle_indices <- sample(nrow(x))
    x_shuffled <- x[shuffle_indices,]
    y_shuffled <- y[shuffle_indices]

    # Assign samples to folds
    folds <- cut(seq(1, nrow(x)), breaks=k, labels=FALSE)

    # Keep track of MSE of every lambda in every fold (matrix of dimension (number of lambdas X number of folds))
    mse_mat <- matrix(0, length(lambdas), k)

    #  Loop through each lambda
    for (i in 1:length(lambdas)) {
        #  Value of lambda:
        i_lambda <- lambdas[i]

        #  Loop through each fold
        for (fold_number in 1:k) {
            #Split data into training data and validation data (based on fold)
            x_train <- x_shuffled[folds != fold_number,]
            y_train <- y_shuffled[folds != fold_number]

            x_test <- x_shuffled[folds == fold_number,]
            y_test <- y_shuffled[folds == fold_number]

            #  Fit the current model using training data
            beta <- fit_model(x_train, y_train, i_lambda, model_type)
            #  Use fitted model to make predictions on testing data
            y_predicted <- model_prediction(x_test, model_type, beta)

            #  Summarize accuracy
            mse <- mean((y_test - y_predicted)^2)
            mse_mat[i, fold_number] <- mse
        }
    }

    #  Compute average(mse) across all folds (lambda_errors is of length(lambdas))
    lambda_errors <- rowMeans(mse_mat)

    #  Compute optimal lambda (lambda with smallest error)
    optimal_lambda <- lambdas[which.min(lambda_errors)]
    return(optimal_lambda)
}


#  Driver to train the glm.
#  x is the expr matrix of dim (samples X genes)
#  y is the corresponding vector of time step of dim (samples)
#  We want to first learn lambda via K-fold CV
#  Then, using the optimal lambda, we will optimize beta
train_glm <- function(x, y, model_type) {

    # Grid of optional lambdas
    lambdas <- c(0, .001, .01, .1, 1)
    
    # Number of folds in CV
    k <- 3
    
    # Search grid of lambdas with k-fold cross validation for optimal lambda
    optimal_lambda <- select_lambda_through_k_fold_cv(x, y, lambdas, k, model_type)

    #  Train model on ALL of training data using optimal lambda
    beta <- fit_model(x, y, optimal_lambda, model_type)
    return(beta)
}





# Compute variance of each row of a matrix (x)
row_var <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

# Return the n top indices of a vector x
get_top_indices <- function(x, n) {
  nx <- length(x)
  p <- nx-n
  xp <- sort(x, partial=p)[p]
  which(x > xp)
}


# Extract indices of the top num_genes with the largest variance.
select_genes_with_largest_variance <- function(rpkm_expr, num_genes) {
    variances <- row_var(rpkm_expr)

    indices <- get_top_indices(variances, num_genes)
    return(indices)
}

#  Visualize our time predictions in the form of a heatmap
plot_heatmap_of_time_predictions <- function(mat, file_name, model_type) {
    #  Convert from matrix form to data frame format
    melted_mat <- melt(mat)
    colnames(melted_mat) <- c("cellLine", "Time","predictedTime")

    #  Use factors to represent cell line and time step
    melted_mat$cellLine <- factor(melted_mat$cellLine)
    melted_mat$trueTime <- factor(melted_mat$Time)

    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=cellLine, y=Time)) + geom_tile(aes(fill=predictedTime)) + scale_fill_gradient2(midpoint=8, guide="colorbar")
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5))
    heatmap <- heatmap + ggtitle(paste0("    ",model_type))

    ggsave(heatmap, file=file_name,width = 15,height=13.5,units="cm")

}

time_step_prediction <- function(sample_info, quant_expr, rpkm_expr, model_type, num_genes) {
    # Extract indices of the top num_genes with the largest variance.
    gene_indices <- select_genes_with_largest_variance(rpkm_expr, num_genes)

    #  Filter expression data to only include those top num_genes genes
    quant_expr_filt <- t(quant_expr[gene_indices,])
    rpkm_expr_filt <- t(rpkm_expr[gene_indices,])

    #  Extract unique cell lines
    cell_lines <- sort(unique(sample_info$cell_line))

    #  Extract unique time steps
    time_steps <- sort(unique(sample_info$time))

    # Keep track of time predictions for each (cell_line,time_step)
    time_prediction_mat <- matrix(NA, length(cell_lines), length(time_steps))
    rownames(time_prediction_mat) <- cell_lines
    colnames(time_prediction_mat) <- time_steps

    #  Loop through each cell line
    for (i in 1:length(cell_lines)) {
        #  ith cell line
        i_cell_line <- cell_lines[i]

        #  Extract indices of samples used for training
        training_indices <- sample_info$cell_line != i_cell_line
        #  Extract indices of samples used for testing
        testing_indices <- sample_info$cell_line == i_cell_line

        #  Extract training data (expr data and corresponding time step)
        training_expr <- quant_expr_filt[training_indices,]
        training_time <- sample_info$time[training_indices]
        #  Extract testing data (expr data and corresponding time step)
        testing_expr <- quant_expr_filt[testing_indices,]
        testing_time <- sample_info$time[testing_indices]

        #  Train glm. This include learning:
        #    1. Coefficient vector (beta) through mle (of dimension (num_genes +1,1))
        #    2. Lasso regression parameter (lambda) through cross validation
        beta <- train_glm(training_expr, training_time, model_type)

        #  Use glm to make predictions
        predicted_times <- model_prediction(testing_expr, model_type, beta)

        #  Add time predictions to matrix (to keep track of everything)
        for (t_index in 1:length(testing_time)) {
            actual_time <- testing_time[t_index]
            predicted_time <- predicted_times[t_index]

            time_prediction_mat[i, actual_time + 1] <- predicted_time
        }
    }
    return(time_prediction_mat)
}




#  Main driver for time step prediction
time_step_prediction_heatmap_plot <- function(sample_info, quant_expr, rpkm_expr, model_type, num_genes, output_file) {
    # Make predictions
    time_prediction_mat <- time_step_prediction(sample_info, quant_expr, rpkm_expr, model_type, num_genes)
    #  Visualize our time predictions in the form of a heatmap
    plot_heatmap_of_time_predictions(time_prediction_mat, output_file, model_type)
}


compare_variance_of_actual_vs_predicted_time <- function(sample_info, quant_expr, rpkm_expr, model_type, num_genes, variance_metric, output_file) {
    #  Extract unique cell lines
    cell_lines <- sort(unique(sample_info$cell_line))
    #  Extract unique time steps
    time_steps <- sort(unique(sample_info$time))
    # Make predictions
    time_prediction_mat <- time_step_prediction(sample_info, quant_expr, rpkm_expr, model_type, num_genes)
    
    # Update sample_info matrix to include time predictions
    sample_info_predicted <- sample_info
    #  Loop through all samples
    for (sample_num in 1:dim(sample_info)[1]) {
        # Find row index of time_prediction_mat corresponding to this sample
        t_cell_line_index <- which(sample_info$cell_line[sample_num] == cell_lines)
        # Find column index of time_prediction_mat corresponding to this sample
        t_time_index <- which(sample_info$time[sample_num] == time_steps)
        # Extract predicted time corresponding to this sample
        sample_info_predicted$time[sample_num] <- round(time_prediction_mat[t_cell_line_index, t_time_index])
        if (sample_info_predicted$time[sample_num] > 15) {
            sample_info_predicted$time[sample_num] <- 15
        }
        if (sample_info_predicted$time[sample_num] < 0) {
            sample_info_predicted$time[sample_num] <- 0
        }
    }

    rpkm_expr <- filter_expr_to_non_zero_genes(rpkm_expr, sample_info)
    rpkm_expr <- filter_expr_to_non_zero_genes(rpkm_expr,sample_info_predicted)

    # Initialize vector of distances 
    distances <- c()
    # Initialize vector of time steps
    t_steps <- c()
    # Initialize predicted_vs_actual
    time_type <- c()

    # Need to compute  sqare distance from mean for each time step seperately..
    # so Loop through each time step
    for (t_step in 0:15) {
        #  Get indices of all samples that are from t_step
        t_sample_indices <- sample_info$time == t_step

        #  Subset rpkm_expr matrix so it only contains samples from t_step
        t_rpkm_expr <- rpkm_expr[,t_sample_indices]

        #  compute square distance from the mean for all (gene, sample) pairs from t_step
        if (variance_metric == "square_distance_from_mean") {
            t_disty <- compute_square_distance_from_mean(t_rpkm_expr)
        }
        if (variance_metric == "sdev") {
            t_disty <- compute_sdev(t_rpkm_expr)
        }
        if (variance_metric == "log_square_distance_from_mean") {
            t_disty <- log(compute_square_distance_from_mean(t_rpkm_expr))
        }
        if (variance_metric == "log_sdev") {
            t_disty <- log(compute_sdev(t_rpkm_expr))
        }
        if (variance_metric == "avg_square_distance_from_mean") {
            t_disty <- compute_avg_square_distance_from_mean(t_rpkm_expr)
        }

        distances <- c(distances, t_disty)
        t_steps <- c(t_steps, numeric(length(t_disty)) + t_step)
        time_type <- c(time_type,rep("actual", length(t_disty)))
    }
    for (t_step in 0:15) {
        #  Get indices of all samples that are from t_step
        t_sample_indices <- sample_info_predicted$time == t_step

        #  Subset rpkm_expr matrix so it only contains samples from t_step
        t_rpkm_expr <- rpkm_expr[,t_sample_indices]

        #  compute square distance from the mean for all (gene, sample) pairs from t_step
        if (variance_metric == "square_distance_from_mean") {
            t_disty <- compute_square_distance_from_mean(t_rpkm_expr)
        }
        if (variance_metric == "sdev") {
            t_disty <- compute_sdev(t_rpkm_expr)
        }
        if (variance_metric == "log_square_distance_from_mean") {
            t_disty <- log(compute_square_distance_from_mean(t_rpkm_expr))
        }
        if (variance_metric == "log_sdev") {
            t_disty <- log(compute_sdev(t_rpkm_expr))
        }
        if (variance_metric == "avg_square_distance_from_mean") {
            t_disty <- compute_avg_square_distance_from_mean(t_rpkm_expr)
        }

        distances <- c(distances, t_disty)
        t_steps <- c(t_steps, numeric(length(t_disty)) + t_step)
        time_type <- c(time_type,rep("predicted", length(t_disty)))
    }

    #  Convert data into data frame format
    df <- data.frame(distance=distances, time=factor(t_steps),time_type=factor(time_type))
    
    testy <- wilcox.test(distances ~ time_type,data=df)

    print(variance_metric)
    print(mean(df[df$time_type == "predicted",]$distance))
    print(mean(df[df$time_type == "actual",]$distance))
    print(median(df[df$time_type == "predicted",]$distance))
    print(median(df[df$time_type == "actual",]$distance))

    # PLOT
    boxplot <- ggplot(df, aes(x=time,y=distance,fill=time_type)) + geom_boxplot(alpha=.7) + labs(x = "Time Step", y = variance_metric) + scale_x_discrete(name="Time Step") + scale_fill_brewer(palette = "Accent")
    boxplot <- boxplot + theme(text = element_text(size=18)) + ggtitle(paste0("Wilcoxon pvalue = ",testy$p.value))
    ggsave(boxplot, file=output_file,width = 20,height=10.5,units="cm")

}





#  compute average square distance from the mean for all samples from t_step
#  Returns vector of length number of (sample,gene) pairs belonging to t_step
compute_square_distance_from_mean <- function(t_rpkm_expr) {

    #  Number of genes/loci
    L <- dim(t_rpkm_expr)[1]

    #  Number of samples belonging to time step t
    N <- dim(t_rpkm_expr)[2]

    # Initialize vector corresponding to the average_square distance from the mean (each element is one (sample, gene) pair)
    disty <- numeric(N*L)

    # compute mean of each gene/locus
    loci_means <- rowMeans(t_rpkm_expr)

    # Compute distance for each sample,gene pair
    counter <- 1
    for (n in 1:N) {
        for (l in 1:L) {
            disty[counter] <- ((t_rpkm_expr[l,n] - loci_means[l])**2)/(loci_means[l]**2)
            counter <- counter + 1
        }
    }
    return(disty)
}

#  compute average square distance from the mean for all samples from t_step
#  Returns vector of length number of (sample,gene) pairs belonging to t_step
compute_avg_square_distance_from_mean <- function(t_rpkm_expr) {

    #  Number of genes/loci
    L <- dim(t_rpkm_expr)[1]

    #  Number of samples belonging to time step t
    N <- dim(t_rpkm_expr)[2]

    # Initialize vector corresponding to the average_square distance from the mean (each element is one (sample, gene) pair)
    disty <- numeric(N)

    # compute mean of each gene/locus
    loci_means <- rowMeans(t_rpkm_expr)

    # Compute distance for each sample,gene pair
    counter <- 1
    for (n in 1:N) {
        disty[counter] <- (N/(L*(N-1)))*(sum(((t_rpkm_expr[,n] - loci_means)**2)/(loci_means**2)))
        counter <- counter + 1
    }
    return(disty)
}

compute_sdev <- function(t_rpkm_expr) {

    #  Number of genes/loci
    L <- dim(t_rpkm_expr)[1]

    #  Number of samples belonging to time step t
    N <- dim(t_rpkm_expr)[2]

    # Initialize vector corresponding to the average_square distance from the mean (each element is one (sample, gene) pair)
    disty <- numeric(L)

    # Compute variance for each gene pair
    counter <- 1
    for (l in 1:L) {
        disty[counter] <- sd(t_rpkm_expr[l,])
        counter <- counter + 1
    }
    return(disty)
}


# Filter to only genes that do not have rpkm == 0 for all sample subsets
# A sample subset is all the samples at time step t
filter_expr_to_non_zero_genes <- function(rpkm_expr_all, sample_info) {
    #  Number of samples
    N <- dim(rpkm_expr_all)[2]
    #  Number of genes
    D <- dim(rpkm_expr_all)[1]

    #  Initialize TRUE/FALSE vector to all TRUE
    #  TRUE represents sample is non-zero expression in each expr matrix subset
    binary_vec <- rep(TRUE, D)

    #  Need to loop through each expr matrix subset (so loop through time steps)
    for (t_step in 0:15) {
        # Sample indices of tth time step
        t_sample_indices <- sample_info$time == t_step

        #  Subseted expression matrix for t time step
        t_rpkm_subset <- rpkm_expr_all[,t_sample_indices]

        # Determine which genes have mapped reads for this sample subset
        t_binary_vec <- !apply(t_rpkm_subset, 1, function(x){all(x==0)})

        binary_vec <- binary_vec*t_binary_vec
    }
    boolean_vec <- binary_vec == 1.0
    return(rpkm_expr_all[boolean_vec,])
}


# Compute the  square distance from mean across all samples (in each time step seperately). Then plot distribution for each time
square_distance_from_mean_driver <- function(sample_info, rpkm_expr_all, metric, square_distance_from_mean_output_file) {
    # Filter to only genes that do not have rpkm == 0 for all sample subsets
    # A sample subset is all the samples at time step t
    rpkm_expr <- filter_expr_to_non_zero_genes(rpkm_expr_all, sample_info)

    # Initialize vector of distances 
    distances <- c()
    # Initialize vector of time steps
    t_steps <- c()

    # Need to compute  sqare distance from mean for each time step seperately..
    # so Loop through each time step
    for (t_step in 0:15) {
        #  Get indices of all samples that are from t_step
        t_sample_indices <- sample_info$time == t_step

        #  Subset rpkm_expr matrix so it only contains samples from t_step
        t_rpkm_expr <- rpkm_expr[,t_sample_indices]

        #  compute square distance from the mean for all (gene, sample) pairs from t_step
        if (metric == "square_distance_from_mean") {
            t_disty <- compute_square_distance_from_mean(t_rpkm_expr)
        }
        if (metric == "sdev") {
            t_disty <- compute_sdev(t_rpkm_expr)
        }
        if (metric == "log_square_distance_from_mean") {
            t_disty <- log(compute_square_distance_from_mean(t_rpkm_expr))
        }
        if (metric == "log_sdev") {
            t_disty <- log(compute_sdev(t_rpkm_expr))
        }
        if (metric == "avg_square_distance_from_mean") {
            t_disty <- compute_avg_square_distance_from_mean(t_rpkm_expr)
        }

        distances <- c(distances, t_disty)
        t_steps <- c(t_steps, numeric(length(t_disty)) + t_step)
    }

    #  Convert data into data frame format
    df <- data.frame(distance=distances, time=factor(t_steps))

    # Plot!

    violin <- ggplot(df, aes(x=time,y=distance)) + geom_boxplot() + labs(x = "Time Step", y = metric)
    violin <- violin + theme(text = element_text(size=18))
    ggsave(violin, file=square_distance_from_mean_output_file,width = 15,height=10.5,units="cm")

}



covariate_pc_pve_heatmap <- function(pc_file, covariate_file, output_file) {
    # Load in data
    pcs <- read.table(pc_file, header=TRUE)
    covs <- read.table(covariate_file, header=TRUE)

    # Remove unimportant columns
    pcs <- as.matrix(pcs[,2:dim(pcs)[2]])
    covs <- data.frame(as.matrix(covs[,3:dim(covs)[2]]))
    # Get covariates into propper class (only necessary for numeric)
    covs$time <- as.numeric(as.character(covs$time))
    covs$library_size <- as.numeric(as.character(covs$library_size))
    covs$feeder_passage <- as.numeric(as.character(covs$feeder_passage))
    covs$feeder_free_passage <- as.numeric(as.character(covs$feeder_free_passage))
    covs$rna_extraction_conc <- as.numeric(as.character(covs$rna_extraction_conc))
    covs$RIN <- as.numeric(as.character(covs$RIN))
    covs$volume_needed <- as.numeric(as.character(covs$volume_needed))
    covs$water_added <- as.numeric(as.character(covs$water_added))
    covs$beating <- as.numeric(as.character(covs$beating))
    covs$line_beating <- as.numeric(as.character(covs$line_beating))
    covs$flash_freezing <- as.numeric(as.character(covs$flash_freezing))
    covs$fastq_percent_duplicates <- as.numeric(as.character(covs$fastq_percent_duplicates))
    covs$percent_gc <- as.numeric(as.character(covs$percent_gc))
    covs$average_sequence_length <- as.numeric(as.character(covs$average_sequence_length))
    covs$total_sequences <- as.numeric(as.character(covs$total_sequences))
    covs$percent_fails <- as.numeric(as.character(covs$percent_fails))


    # Initialize PVE heatmap
    pve_map <- matrix(0, dim(covs)[2], dim(pcs)[2])
    colnames(pve_map) <- colnames(pcs)
    rownames(pve_map) <- colnames(covs)

    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(pcs)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_cov in 1:num_covs) {
            pc_vec <- pcs[,num_pc]
            cov_vec <- covs[,num_cov]
            lin_model <- lm(pc_vec ~ cov_vec)
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    
    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Covariate", "PC","PVE")

    #  Use factors to represent covariate and pc name
    melted_mat$Covariate <- factor(melted_mat$Covariate, levels = rownames(pve_map)[ord])
    melted_mat$PC <- factor(melted_mat$PC)
    levels(melted_mat$PC) <- c(levels(melted_mat$PC)[1],levels(melted_mat$PC)[3:10],levels(melted_mat$PC)[2])



    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Covariate, y=PC)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + theme(text = element_text(size=14), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5))

    # Save File
    ggsave(heatmap, file=output_file,width = 19,height=13.5,units="cm")
}

covariate_pc_specific_genes_pve_heatmap <- function(pc_file, quant_expr, covariate_file, output_file) {
    # Load in data
    pcs <- read.table(pc_file, header=TRUE)
    covs <- read.table(covariate_file, header=TRUE)

    # Remove unimportant columns
    pcs <- as.matrix(pcs[,2:dim(pcs)[2]])
    covs <- data.frame(as.matrix(covs[,3:dim(covs)[2]]))
    # Get covariates into propper class (only necessary for numeric)
    covs$time <- as.numeric(as.character(covs$time))
    covs$library_size <- as.numeric(as.character(covs$library_size))
    covs$feeder_passage <- as.numeric(as.character(covs$feeder_passage))
    covs$feeder_free_passage <- as.numeric(as.character(covs$feeder_free_passage))
    covs$rna_extraction_conc <- as.numeric(as.character(covs$rna_extraction_conc))
    covs$RIN <- as.numeric(as.character(covs$RIN))
    covs$volume_needed <- as.numeric(as.character(covs$volume_needed))
    covs$water_added <- as.numeric(as.character(covs$water_added))
    covs$beating <- as.numeric(as.character(covs$beating))
    covs$line_beating <- as.numeric(as.character(covs$line_beating))
    covs$flash_freezing <- as.numeric(as.character(covs$flash_freezing))
    covs$fastq_percent_duplicates <- as.numeric(as.character(covs$fastq_percent_duplicates))
    covs$percent_gc <- as.numeric(as.character(covs$percent_gc))
    covs$average_sequence_length <- as.numeric(as.character(covs$average_sequence_length))
    covs$total_sequences <- as.numeric(as.character(covs$total_sequences))
    covs$percent_fails <- as.numeric(as.character(covs$percent_fails))

    # Add specific gene info
    row_label <- which(rownames(quant_expr) == "ENSG00000118194")
    quant_expr <- as.matrix(quant_expr)
    covs$Troponin <- as.vector(quant_expr[row_label,])

    row_label <- which(rownames(quant_expr) == "ENSG00000181449")
    covs$Sox2 <- as.vector(quant_expr[row_label,])

    # Initialize PVE heatmap
    pve_map <- matrix(0, dim(covs)[2], dim(pcs)[2])
    colnames(pve_map) <- colnames(pcs)
    rownames(pve_map) <- colnames(covs)

    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(pcs)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_cov in 1:num_covs) {
            pc_vec <- pcs[,num_pc]
            cov_vec <- covs[,num_cov]
            lin_model <- lm(pc_vec ~ cov_vec)
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    
    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Covariate", "PC","PVE")

    #  Use factors to represent covariate and pc name
    melted_mat$Covariate <- factor(melted_mat$Covariate, levels = rownames(pve_map)[ord])
    melted_mat$PC <- factor(melted_mat$PC)
    levels(melted_mat$PC) <- c(levels(melted_mat$PC)[1],levels(melted_mat$PC)[3:10],levels(melted_mat$PC)[2])



    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Covariate, y=PC)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + theme(text = element_text(size=14), panel.background = element_blank(), axis.text.x = element_text(angle = 90, vjust=.5))

    # Save File
    ggsave(heatmap, file=output_file,width = 19,height=13.5,units="cm")
}



#####################################################################################################
# Load Data
#####################################################################################################

options(warn=1)
preprocess_total_expression_dir = args[1]  # Where total expression processed data is_autosomal
visualize_total_expression_dir = args[2]  # Ouputdir to save images
covariate_dir = args[3]  # Input dir with covariate information

#  Get sample information 
sample_info_file <- paste0(preprocess_total_expression_dir, "sample_info.txt")
sample_info <- read.table(sample_info_file, header=TRUE)

#  Get quantile normalized expression data
quantile_normalized_exp_file <- paste0(preprocess_total_expression_dir, "quantile_normalized.txt")
quant_expr <- read.csv(quantile_normalized_exp_file, header=TRUE, sep=" ")

#  Get rpkm expression_data
rpkm_exp_file <- paste0(preprocess_total_expression_dir, "rpkm.txt")
rpkm_expr <- read.csv(rpkm_exp_file, header=TRUE, sep=" ")

#  Get covariate file
covariate_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
covariates <- read.table(covariate_file,header=TRUE)



#####################################################################################################
# Run Analysis / Create Plots
#####################################################################################################



####################################################################
# Plot library size
####################################################################

################
# Make barplot showing library sizes of each sample
library_size_output_file <- paste0(visualize_total_expression_dir, "library_size.pdf")
plot_library_size(sample_info, library_size_output_file)



####################################################################
# PCA plots with samples labeled by various covariates
####################################################################

##################
#  Perform PCA. Plot first 2 pcs as a function of time step 
pca_plot_time_step_output_file <- paste0(visualize_total_expression_dir, "pca_plot_1_2_time_step.pdf")
plot_pca_time_step(sample_info, quant_expr, pca_plot_time_step_output_file)



#################
# Perform PCA. Plot specified PCs as a function of gene expression of:
# a. Troponin (gene expressed in cardiomyocytes)
# b. sox2 (gene expressed in ipscs)
pc_num1 <- 1
pc_num2 <- 2
ensamble_id <- "ENSG00000118194"
gene_name <- "Troponin"
pca_plot_gene_filled_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_",gene_name,"_gene_filled.pdf")
plot_pca_real_valued_gene_filled(sample_info, quant_expr, ensamble_id,gene_name,pc_num1,pc_num2,pca_plot_gene_filled_output_file)

pc_num1 <- 2
pc_num2 <- 3
ensamble_id <- "ENSG00000118194"
gene_name <- "Troponin"
pca_plot_gene_filled_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_",gene_name,"_gene_filled.pdf")
plot_pca_real_valued_gene_filled(sample_info, quant_expr, ensamble_id,gene_name,pc_num1,pc_num2,pca_plot_gene_filled_output_file)

pc_num1 <- 1
pc_num2 <- 2
ensamble_id <- "ENSG00000181449"
gene_name <- "sox2"
pca_plot_gene_filled_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_",gene_name,"_gene_filled.pdf")
plot_pca_real_valued_gene_filled(sample_info, quant_expr, ensamble_id,gene_name,pc_num1,pc_num2,pca_plot_gene_filled_output_file)

pc_num1 <- 2
pc_num2 <- 3
ensamble_id <- "ENSG00000181449"
gene_name <- "sox2"
pca_plot_gene_filled_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_",gene_name,"_gene_filled.pdf")
plot_pca_real_valued_gene_filled(sample_info, quant_expr, ensamble_id,gene_name,pc_num1,pc_num2,pca_plot_gene_filled_output_file)



#################
# Perform PCA. Plot specified PCS as a function of various covariates:
# a. cell line
# b. rna extraction person
# c. rna extraction round
# d. differentiation batch

pc_num1 <- 1
pc_num2 <- 2

pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_cell_line.pdf")
plot_pca_categorical_covariate(sample_info, quant_expr, pca_plot_cell_line_output_file,factor(sample_info$cell_line), "cell_line", pc_num1,pc_num2)

pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_rna_extraction_persion.pdf")
plot_pca_categorical_covariate(sample_info, quant_expr, pca_plot_cell_line_output_file,factor(covariates$RNA_extraction_person), "rna_extraction_person", pc_num1,pc_num2)

pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_rna_extraction_round.pdf")
plot_pca_categorical_covariate(sample_info, quant_expr, pca_plot_cell_line_output_file,factor(covariates$RNA_extraction_round), "rna_extraction_round", pc_num1,pc_num2)

pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_differentiation_batch.pdf")
plot_pca_categorical_covariate(sample_info, quant_expr, pca_plot_cell_line_output_file,factor(covariates$differentiation_batch), "differentiation_batch", pc_num1,pc_num2)


#################
# Perform PCA. Make seperate plot for each cell line:

pc_num1<-1
pc_num2<-2
pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_seperate_cell_lines.pdf")
plot_pca_seperate_cell_lines(sample_info, quant_expr, pca_plot_cell_line_output_file,pc_num1,pc_num2)

pc_num1<-2
pc_num2<-3
pca_plot_cell_line_output_file <- paste0(visualize_total_expression_dir, "pca_plot_",pc_num1,"_",pc_num2,"_seperate_cell_lines.pdf")
plot_pca_seperate_cell_lines(sample_info, quant_expr, pca_plot_cell_line_output_file,pc_num1,pc_num2)

#################
# Perform PCA. Plot variance explained of the first n PCs:
n <- 20
pca_plot_variance_explained_output_file <- paste0(visualize_total_expression_dir, "pca_plot_variance_explained", n, ".pdf")
plot_pca_variance_explained(sample_info, quant_expr, n, pca_plot_variance_explained_output_file)





####################################################################
# Time course of expression of specific genes seperated by cell line
####################################################################
ensamble_id <- "ENSG00000118194"
gene_name <- "Troponin"
line_plot_file <- paste0(visualize_total_expression_dir, gene_name,"_time_course_grouped_by_cell_line.pdf")
gene_time_course_line_plot_grouped_by_cell_line(sample_info, quant_expr, ensamble_id, gene_name, line_plot_file)

ensamble_id <- "ENSG00000181449"
gene_name <- "sox2"
line_plot_file <- paste0(visualize_total_expression_dir, gene_name,"_time_course_grouped_by_cell_line.pdf")
gene_time_course_line_plot_grouped_by_cell_line(sample_info, quant_expr, ensamble_id, gene_name, line_plot_file)




###########################################################################################
# Scatter plot of 2nd (and others) PC loading vs a genes expression levels in t = 15 samples
###########################################################################################
ensamble_id <- "ENSG00000118194"
gene_name <- "Troponin"
time_step <- 15
pc_num <- 2

pc_gene_scatter_output_file <- paste0(visualize_total_expression_dir, "pc_num_", pc_num, "_",gene_name,"_time_step_",time_step,"_scatter.pdf")
pc_gene_scatter(sample_info, quant_expr, ensamble_id, gene_name, time_step, pc_num, pc_gene_scatter_output_file)


ensamble_id <- "ENSG00000181449"
gene_name <- "sox2"
time_step <- 15
pc_num <- 2

pc_gene_scatter_output_file <- paste0(visualize_total_expression_dir, "pc_num_", pc_num, "_",gene_name,"_time_step_",time_step,"_scatter.pdf")
pc_gene_scatter(sample_info, quant_expr, ensamble_id, gene_name, time_step, pc_num, pc_gene_scatter_output_file)




####################################################################
# Covariates explaining variance in principle components
####################################################################


# Make heatmap showing PVE between pcs and covariates
pc_file <- paste0(covariate_dir,"principal_components_10.txt")
covariate_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
output_file <- paste0(visualize_total_expression_dir, "pc_covariate_pve_heatmap.png")
covariate_pc_pve_heatmap(pc_file, covariate_file,output_file)

# Make heatmap showing PVE between pcs & (covariates and troponin/sox2 expression)
pc_file <- paste0(covariate_dir,"principal_components_10.txt")
covariate_file <- paste0(covariate_dir, "processed_covariates_categorical.txt")
output_file <- paste0(visualize_total_expression_dir, "pc_covariate_troponin_sox2_pve_heatmap.png")
covariate_pc_specific_genes_pve_heatmap(pc_file, quant_expr, covariate_file,output_file)



####################################################################
# Time step Predictions (ignore for now)
####################################################################
#  Predict time step of each sample using sparse linear regression model
model_type <- "linear_regression"  # Options are "linear_regression" and "beta_binomial_regression"
num_genes <- 1000  #  Filter genes to top num_genes with largest variance
output_file <- paste0(visualize_total_expression_dir, "time_prediction_", model_type, "_", num_genes, ".pdf")
time_step_prediction_heatmap_plot(sample_info, quant_expr, rpkm_expr, model_type, num_genes, output_file)


#  Predict time step of each sample using sparse beta-binomial generalized linear model
model_type <- "beta_binomial_regression"  # Options are "linear_regression" and "beta_binomial_regression"
num_genes <- 1000  #  Filter genes to top num_genes with largest variance
output_file <- paste0(visualize_total_expression_dir, "time_prediction_", model_type, "_", num_genes, ".pdf")
time_step_prediction_heatmap_plot(sample_info, quant_expr, rpkm_expr, model_type, num_genes, output_file)






####################################################################
# Various metrics to access the variance of the samples (binned by both time step and actual time vs pseudo time)
####################################################################

# Compute the square distance from mean across all (gene,samples) (in each time step seperately). Then plot distribution for each time
metric<-"log_square_distance_from_mean"
square_distance_from_mean_output_file <- paste0(visualize_total_expression_dir, metric, ".pdf")
#square_distance_from_mean_driver(sample_info, rpkm_expr, metric, square_distance_from_mean_output_file)
# Compute the sdev of each gene (in each time step seperately). Then plot distribution for each time
metric<-"log_sdev"
square_distance_from_mean_output_file <- paste0(visualize_total_expression_dir, metric, ".pdf")
#square_distance_from_mean_driver(sample_info, rpkm_expr, metric, square_distance_from_mean_output_file)
# Compute the average square distance from mean across all samples(in each time step seperately). Then plot distribution for each time
metric<-"avg_square_distance_from_mean"
square_distance_from_mean_output_file <- paste0(visualize_total_expression_dir, metric, ".pdf")
#square_distance_from_mean_driver(sample_info, rpkm_expr, metric, square_distance_from_mean_output_file)


#  Compute predicted time steps of each sample using sparse linear regression model.
#  Then compare variance under predicted time steps vs actual time steps
model_type <- "linear_regression" # Options are "linear_regression" and "beta_binomial_regression"
num_genes <- 1000 # Filter genes to top num_genes with largest variance
variance_metric <- "log_sdev"
output_file <- paste0(visualize_total_expression_dir, "actual_vs_predicted_time_",model_type,"_",num_genes,"_", variance_metric, ".pdf")
#compare_variance_of_actual_vs_predicted_time(sample_info, quant_expr, rpkm_expr, model_type, num_genes, variance_metric, output_file)

model_type <- "linear_regression" # Options are "linear_regression" and "beta_binomial_regression"
num_genes <- 1000 # Filter genes to top num_genes with largest variance
variance_metric <- "avg_square_distance_from_mean"
output_file <- paste0(visualize_total_expression_dir, "actual_vs_predicted_time_",model_type,"_",num_genes,"_", variance_metric, ".pdf")
#compare_variance_of_actual_vs_predicted_time(sample_info, quant_expr, rpkm_expr, model_type, num_genes, variance_metric, output_file)

model_type <- "linear_regression" # Options are "linear_regression" and "beta_binomial_regression"
num_genes <- 1000 # Filter genes to top num_genes with largest variance
variance_metric <- "log_square_distance_from_mean"
output_file <- paste0(visualize_total_expression_dir, "actual_vs_predicted_time_",model_type,"_",num_genes,"_", variance_metric, ".pdf")
#compare_variance_of_actual_vs_predicted_time(sample_info, quant_expr, rpkm_expr, model_type, num_genes, variance_metric, output_file)
#  Compute predicted time steps of each sample using sparse linear regression model.
#  Then compare variance under predicted time steps vs actual time steps
model_type <- "beta_binomial_regression" # Options are "linear_regression" and "beta_binomial_regression"
num_genes <- 1000 # Filter genes to top num_genes with largest variance
variance_metric <- "log_sdev"
output_file <- paste0(visualize_total_expression_dir, "actual_vs_predicted_time_",model_type,"_",num_genes,"_", variance_metric, ".pdf")
#compare_variance_of_actual_vs_predicted_time(sample_info, quant_expr, rpkm_expr, model_type, num_genes, variance_metric, output_file)

model_type <- "beta_binomial_regression" # Options are "linear_regression" and "beta_binomial_regression"
num_genes <- 1000 # Filter genes to top num_genes with largest variance
variance_metric <- "avg_square_distance_from_mean"
output_file <- paste0(visualize_total_expression_dir, "actual_vs_predicted_time_",model_type,"_",num_genes,"_", variance_metric, ".pdf")
#compare_variance_of_actual_vs_predicted_time(sample_info, quant_expr, rpkm_expr, model_type, num_genes, variance_metric, output_file)

model_type <- "beta_binomial_regression" # Options are "linear_regression" and "beta_binomial_regression"
num_genes <- 1000 # Filter genes to top num_genes with largest variance
variance_metric <- "log_square_distance_from_mean"
output_file <- paste0(visualize_total_expression_dir, "actual_vs_predicted_time_",model_type,"_",num_genes,"_", variance_metric, ".pdf")
#compare_variance_of_actual_vs_predicted_time(sample_info, quant_expr, rpkm_expr, model_type, num_genes, variance_metric, output_file)
