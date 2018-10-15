args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(reshape)
library(cowplot)





# Plot correlation histogram for summary stat
symmetric_correlation_heatmap_general <- function(correlation_matrix, output_file) {
    nn <- dim(correlation_matrix)[1]
    vec <- c()
    for (i in 1:nn){
        for (j in 1:nn) {
            if (i != j) {
                vec <- c(vec,correlation_matrix[i,j])
            }
        }
    }
    maxy <- max(vec)
    for (i in 1:nn) {
        correlation_matrix[i,i] <- maxy
    }
    melted_corr <- melt(correlation_matrix)

    ord <- hclust( dist(scale(correlation_matrix), method = "euclidean"), method = "ward.D" )$order

    print(ord)


    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1, levels = rownames(correlation_matrix)[ord])
    melted_corr$X2 <- factor(melted_corr$X2, levels = colnames(correlation_matrix)[ord])




    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    #heatmap <- heatmap + scale_fill_distiller()
    #heatmap <- heatmap + scale_fill_brewer(values = brewer.pal(3,"RdPu"))
    heatmap <- heatmap + scale_fill_distiller(palette = "RdPu", direction=1)
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 90))
    heatmap <- heatmap + labs(x = "Sample ID", y = "Sample ID", fill= "Spearman Rho")
    heatmap <- heatmap +   theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

    ggsave(heatmap, file=output_file,width = 50,height=36,units="cm")


}







output_dir = args[1] 

input_file <- paste0(output_dir, "raw_counts_cmp_data_sets.txt")
output_file <- paste0(output_dir, "raw_counts_cmp_data_sets_heatmap.png")

raw_data <- read.table(input_file,header=TRUE)
counts <- raw_data[,2:(dim(raw_data)[2])]
corr_mat <- cor(counts,method="spearman")

symmetric_correlation_heatmap_general(corr_mat, output_file)
