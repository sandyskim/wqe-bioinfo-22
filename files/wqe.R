library('data.table')
library('nimble')
library('dplyr')
library('biomaRt')
library('rstan')
library('ggplot2')
library('cowplot')
setwd("~/Desktop")

gtex_gene_counts <- fread(file="/Users/sandykim/Desktop/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'https://www.ensembl.org')
genes <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22)), mart = mart)

gtex_gene_counts$Name <- gsub("\\..*","", gtex_gene_counts$Name)
gtex_gene_counts <- gtex_gene_counts[which(gtex_gene_counts$Name %in% genes$ensembl_gene_id),]
gtex_gene_counts[which(gtex_gene_counts$Name == 'ENSG00000283992'),]$Description = 'SLURP2'

gtex_gene_counts <- gtex_gene_counts[,3:17384]
gtex_gene_counts <- gtex_gene_counts[1:100, 1:100]

gtex_gene_counts <- gtex_gene_counts[which(apply(gtex_gene_counts, 1, var) != 0),]
pca_counts <- prcomp(t(gtex_gene_counts), center=TRUE, scale=TRUE)
pca_counts <- as.data.frame(pca_counts$x)
pca_counts <- (data.frame(pca_counts[,1:2], outlier=pca_counts$PC1 < -8))
p <- ggplot(data=pca_counts) + geom_point(aes(x=PC1, y=PC2, color=outlier), alpha=0.5) + labs(title = "PCA biplot, GTEx read counts of 100 samples, 100 genes", color="outlier")
p

which(pca_counts$outlier == TRUE)
save_plot(p, file="pca_gtex.png", base_width = 14, base_height = 10)

mus <- fread("gtex_mus.tsv")


set.seed(42)
dispersion <- 100
n_samples <- 100
n_genes <- 100
counts <- as.data.frame(matrix(nrow=n_genes, ncol=n_samples))
p_outliers <- 1.0
n_outliers <- round(n_samples * p_outliers)
outlier <- rep(0, times=n_samples)
outlier[1:n_outliers] <- sample(c(-1, 1), n_outliers, replace = TRUE)
outlier_effect_p_1.0 <- exp(outlier)
for (i in (1:n_samples)) {
  counts[,i] <- rnbinom(n_genes, size = dispersion, mu = mus$mu * outlier_effect[i])
}

outlier_effects <- data.frame("0.0" = outlier_effect,
                              "0.1" = outlier_effect_p_0.1,
                              "0.2" = outlier_effect_p_0.2,
                              "0.3" = outlier_effect_p_0.3,
                              "0.4" = outlier_effect_p_0.4,
                              "0.5" = outlier_effect_p_0.5,
                              "0.6" = outlier_effect_p_0.6,
                              "0.7" = outlier_effect_p_0.7,
                              "0.8" = outlier_effect_p_0.8,
                              "0.9" = outlier_effect_p_0.9,
                              "1.0" = outlier_effect_p_1.0)

saveRDS(outlier_effects, file="outlier_effects.rds")

rownames(counts) <- mus$V1[1:n_genes]

counts <- counts[which(apply(counts, 1, var) != 0),]
pca_counts <- prcomp(t(counts), center=TRUE, scale=TRUE)

p0 <- ggplot(data=pca_counts_p_0.0) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 0.0", color = "outlier")
p1 <- ggplot(data=pca_counts_p_0.1) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 0.1", color = "outlier")
p2 <- ggplot(data=pca_counts_p_0.2) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 0.2", color = "outlier")
p3 <- ggplot(data=pca_counts_p_0.3) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 0.3", color = "outlier")
p4 <- ggplot(data=pca_counts_p_0.4) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 0.4", color = "outlier")
p5 <- ggplot(data=pca_counts_p_0.5) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 0.5", color = "outlier")
p6 <- ggplot(data=pca_counts_p_0.6) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 0.6", color = "outlier")
p7 <- ggplot(data=pca_counts_p_0.7) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 0.7", color = "outlier")
p8 <- ggplot(data=pca_counts_p_0.8) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 0.8", color = "outlier")
p9 <- ggplot(data=pca_counts_p_0.9) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 0.9", color = "outlier")
p10 <- ggplot(data=pca_counts_p_1.0) + geom_point(aes(x=PC1, y=PC2, color=as.logical(outlier)), alpha=0.5) + labs(title = "p_outlier = 1.0", color = "outlier")

p11 <- plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, vjust = -0.5)
title <- ggdraw() + draw_label("PCA biplot, simulated RNA-seq read counts of 100 samples, 100 genes", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p12 <- plot_grid(title, p11, ncol = 1,  rel_heights = c(.1,.9))
p12

save_plot(p12, file="pca.png", base_width = 14, base_height = 10)


ggsave("simulation.png", p, width=14, height=10)

wqeConsts <- list(
  n_samples = ncol(gtex_gene_counts),
  n_genes = nrow(gtex_gene_counts),
  mean_counts = apply(gtex_gene_counts, 1, mean)
)

wqeData <-  list(counts = gtex_gene_counts)

wqeModel <- nimbleCode({
  pi ~ dbeta(1, 10)
  sigma ~ dgamma(1, 0.1)
  phi ~ dgamma (1, 0.1)
  tau ~ dgamma(1, 1)
  
  for (n in 1:n_samples) {
    psi[n] ~ dbern(pi)
    outlier[n] ~ dnorm(0, sigma)
    log_outlier_effect[n] <- psi[n] * dnorm(outlier[n], tau)
    outlier_effect[n] <- 2^(log_outlier_effect[n])
  }
  for (g in 1:n_genes) {
    phi_genes[g] ~ dgamma(phi, 1);
    for (n in 1:n_samples) {
      mu[g,n] <- mean_counts[g] * outlier_effect[n]
      var[g,n] <- mu[g,n] + (mu[g,n]^2/phi_genes[g]) 
      p[g,n] <- (var[g,n] - mu[g,n])/var[g,n]
      r[g,n] <- (mu[g,n]^2)/(var[g,n] - mu[g,n])
      counts[g,n] ~ dnegbin(p[g,n], r[g,n])
    }
  }
})

wqe <- nimbleModel(code = wqeModel,
                      constants = wqeConsts,
                      data = wqeData)

mcmc_configuration = configureMCMC(wqe, print=TRUE, monitors=c('phi', 'phi_genes', 'pi', 'tau', 'sigma', 'psi', 'outlier_effect', 'log_outlier_effect'))

wqe_mcmc <- buildMCMC(mcmc_configuration)
c_wqe_mcmc <- compileNimble(wqe)
c_wqe_mcmc <- compileNimble(wqe_mcmc, project = wqe)

mcmc_samples <- runMCMC(c_wqe_mcmc,
                        nchains=2,
                        niter=2000,
                        nburnin=1000,
                        summary=TRUE)

mcmc_summary <- as.data.frame(mcmc_samples$summary$all.chains)

mcmc_summary_p_0.0 <- mcmc_summary
mcmc_summaries <- data.frame(mcmc_summary_p_0.1, mcmc_summary_p_0.2, mcmc_summary_p_0.3, mcmc_summary_p_0.4, mcmc_summary_p_0.5,
                    mcmc_summary_p_0.6, mcmc_summary_p_0.7, mcmc_summary_p_0.8, mcmc_summary_p_0.9, mcmc_summary_p_1.0)

for (i in length(mcmc_summaries)) {
  saveRDS(paste0("mcmc_summary_p_", format(round(i, digits=2), nsmall=1)), paste0("mcmc_summary_p_", format(round(i, digits=2), nsmall=1),".rds"))
}
saveRDS(mcmc_summary_p_1.0, "mcmc_summary_p_1.0.rds")
mcmc_summary <- mcmc_summary_p_1.0
mcmc_summary_p_0.0 <- readRDS("mcmc_summary_p_0.0.rds")
mcmc_summary_p_0.2 <- readRDS("mcmc_summary_p_0.2.rds")
mcmc_summary_p_0.4 <- readRDS("mcmc_summary_p_0.4.rds")
mcmc_summary_p_0.6 <- readRDS("mcmc_summary_p_0.6.rds")
mcmc_summary_p_0.8 <- readRDS("mcmc_summary_p_0.8.rds")
mcmc_summary_p_1.0 <- readRDS("mcmc_summary_p_1.0.rds")

oe_gtex = data.frame(mcmc_summary[grepl('outlier_effect', rownames(mcmc_summary)), ][101:200,])
outlier_gtex <- pca_counts$outlier
truth_table_gtex <- as.data.frame(table(reference = outlier_gtex, prediction = as.logical(round(oe_gtex$Mean-1,4))))
p0 <- ggplot(truth_table_gtex, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev) + labs(title = "GTEx confusion matrix") 
p0

save_plot(p0, file="confusion_gtex.png", base_width = 14, base_height = 10)


p1 <- ggplot(data=data.frame(outlier_gtex)) +   geom_point(aes(x=1:100, y=outlier_gtex), alpha=0.5) +
  labs(title="simulated outlier effect sizes", x="sample", y = "outlier effect size") + geom_vline(xintercept=0, alpha=0.25)
p1
p2 <- ggplot(data=data.frame(oe_gtex)) +   geom_point(aes(x=1:100, y=Mean), alpha=0.5) +
  labs(title="simulated outlier effect sizes", x="sample", y = "outlier effect size") + geom_vline(xintercept=0, alpha=0.25)
p2

p3 <- plot_grid(p1, p2, align="v", ncol=1)
p3

oe_0.0 = data.frame(mcmc_summary_p_0.0[grepl('outlier_effect', rownames(mcmc_summary)), ][101:200,])
oe_0.2 = data.frame(mcmc_summary_p_0.2[grepl('outlier_effect', rownames(mcmc_summary)), ][101:200,])
oe_0.4 = data.frame(mcmc_summary_p_0.4[grepl('outlier_effect', rownames(mcmc_summary)), ][101:200,])
oe_0.6 = data.frame(mcmc_summary_p_0.6[grepl('outlier_effect', rownames(mcmc_summary)), ][101:200,])
oe_0.8 = data.frame(mcmc_summary_p_0.8[grepl('outlier_effect', rownames(mcmc_summary)), ][101:200,])
oe_1.0 = data.frame(mcmc_summary_p_1.0[grepl('outlier_effect', rownames(mcmc_summary)), ][101:200,])

p1 <- ggplot(data=data.frame(outlier_effects)) +   geom_point(aes(x=1:100, y=X0.0), alpha=0.5) +
  labs(title="simulated outlier effect sizes", x="sample", y = "outlier effect size") + geom_vline(xintercept=0, alpha=0.25)
p2 <- ggplot(data=round(oe_0.0,2)) + geom_point(aes(x=1:100, y=Mean), alpha=0.5) + 
  labs(title="inferred outlier effect sizes", x="sample", y = "mean outlier effect size") + geom_vline(xintercept=0, alpha=0.25)
title <- ggdraw() + draw_label("p_outlier=0.0", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p3 <- plot_grid(p1, p2, ncol=1, align="v")
p3 <- plot_grid(title, p3, ncol = 1,  rel_heights = c(.1,.9))
p3

p4 <- ggplot(data=data.frame(outlier_effects)) +   geom_point(aes(x=1:100, y=X0.2), alpha=0.5) +
  labs(title="simulated outlier effect sizes", x="sample", y = "outlier effect size") + geom_vline(xintercept=20, alpha=0.25)
p5 <- ggplot(data=oe_0.2) + geom_point(aes(x=1:100, y=Mean), alpha=0.5) + 
  labs(title="inferred outlier effect sizes", x="sample", y = "mean outlier effect size") + geom_vline(xintercept=20, alpha=0.25)
title <- ggdraw() + draw_label("p_outlier=0.2", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p6 <- plot_grid(p4, p5, ncol=1, align="v")
p6 <- plot_grid(title, p6, ncol = 1,  rel_heights = c(.1,.9))
p6

p7 <- ggplot(data=data.frame(outlier_effects)) +   geom_point(aes(x=1:100, y=X0.4), alpha=0.5) +
  labs(title="simulated outlier effect sizes", x="sample", y = "outlier effect size") + geom_vline(xintercept=40, alpha=0.25)
p8 <- ggplot(data=oe_0.4) + geom_point(aes(x=1:100, y=Mean), alpha=0.5) + 
  labs(title="inferred outlier effect sizes", x="sample", y = "mean outlier effect size") + geom_vline(xintercept=40, alpha=0.25)
title <- ggdraw() + draw_label("p_outlier=0.4", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p9 <- plot_grid(p7, p8, ncol=1, align="v")
p9 <- plot_grid(title, p9, ncol = 1,  rel_heights = c(.1,.9))
p9

p10 <- ggplot(data=data.frame(outlier_effects)) +   geom_point(aes(x=1:100, y=X0.6), alpha=0.5) +
  labs(title="simulated outlier effect sizes", x="sample", y = "outlier effect size") + geom_vline(xintercept=60, alpha=0.25)
p11 <- ggplot(data=oe_0.6) + geom_point(aes(x=1:100, y=Mean), alpha=0.5) + 
  labs(title="inferred outlier effect sizes", x="sample", y = "mean outlier effect size") + geom_vline(xintercept=60, alpha=0.25)
title <- ggdraw() + draw_label("p_outlier=0.6", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p12 <- plot_grid(p10, p11, ncol=1, align="v")
p12 <- plot_grid(title, p12, ncol = 1,  rel_heights = c(.1,.9))
p12

p13 <- ggplot(data=data.frame(outlier_effects)) +   geom_point(aes(x=1:100, y=X0.8), alpha=0.5) +
  labs(title="simulated outlier effect sizes", x="sample", y = "outlier effect size") + geom_vline(xintercept=80, alpha=0.25)
p14 <- ggplot(data=oe_0.8) + geom_point(aes(x=1:100, y=Mean), alpha=0.5) + 
  labs(title="inferred outlier effect sizes", x="sample", y = "mean outlier effect size") + geom_vline(xintercept=80, alpha=0.25)
title <- ggdraw() + draw_label("p_outlier=0.8", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p15 <- plot_grid(p13, p14, ncol=1, align="v")
p15 <- plot_grid(title, p15, ncol = 1,  rel_heights = c(.1,.9))
p15

p16 <- ggplot(data=data.frame(outlier_effects)) +   geom_point(aes(x=1:100, y=X1.0), alpha=0.5) +
  labs(title="simulated outlier effect sizes", x="sample", y = "outlier effect size") + geom_vline(xintercept=100, alpha=0.25)
p17 <- ggplot(data=oe_1.0) + geom_point(aes(x=1:100, y=Mean), alpha=0.5) + 
  labs(title="inferred outlier effect sizes", x="sample", y = "mean outlier effect size") + geom_vline(xintercept=100, alpha=0.25)
title <- ggdraw() + draw_label("p_outlier=1.0", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p18 <- plot_grid(p16, p17, ncol=1, align="v")
p18 <- plot_grid(title, p18, ncol = 1,  rel_heights = c(.1,.9))
p18

title <- ggdraw() + draw_label("outlier effect sizes", fontface="bold", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p19 <- plot_grid(p3, p6, p9, p12, p15, p18, ncol=2, align="v")
p19 <- plot_grid(title, p19, ncol = 1,  rel_heights = c(.05,.95))

ggsave("oes.png", p19, width=14, height=10)

psi = mcmc_summary[grepl('psi', rownames(mcmc_summary)), ][101:200,]
psi = data.frame(psi)
plot(psi$Mean)

p <- ggplot(data=outlier_effect)

truth_table <- data.frame(reference = as.logical(outlier))
truth_table <- data.frame(truth_table, prediction = as.logical(round(oe$Mean-1, 2)))
View(oe)

plot(psi$Mean)

alpha = 0.90
n_effects = 0.10 * 100
sensitivity = mean(psi$Mean[1:n_effects] >= 1 - alpha)
sensitivity
fdr = sum(psi$Mean[n_effects:nrow(psi)] >= 1 - alpha) / sum(psi$Mean > 1 - alpha)
fdr
# truth_table <- as.data.frame(table(truth_table))

p0 <- ggplot(truth_table_p_0, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev) + labs(title = "p_outlier = 0.0") 
p1 <- ggplot(truth_table_p_0.1, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev) + labs(title = "p_outlier = 0.1") 
p2 <- ggplot(truth_table_p_0.2, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev)  + labs(title = "p_outlier = 0.2") 
p3 <- ggplot(truth_table_p_0.3, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev)  + labs(title = "p_outlier = 0.3") 
p4 <- ggplot(truth_table_p_0.4, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev)  + labs(title = "p_outlier = 0.4") 
p5 <- ggplot(truth_table_p_0.5, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev) + labs(title = "p_outlier = 0.5") 
p6 <- ggplot(truth_table_p_0.6, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev) + labs(title = "p_outlier = 0.6") 
p7 <- ggplot(truth_table_p_0.7, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev)  + labs(title = "p_outlier = 0.7") 
p8 <- ggplot(truth_table_p_0.8, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev) + labs(title = "p_outlier = 0.8") 
p9 <- ggplot(truth_table_p_0.9, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev)  + labs(title = "p_outlier = 0.9") 
p10 <- ggplot(truth_table_p_1.0, aes(x=prediction, y= reference, fill= Freq)) + geom_tile(alpha=1) + geom_text(aes(label=Freq), colour = "white")  + 
  theme_minimal() + scale_x_discrete(limits = rev) + labs(title = "p_outlier = 1.0") 
library(cowplot)
title <- ggdraw() + draw_label("outlier confusion matrix", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
p11 <- plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
p12 <- plot_grid(title, p11, ncol = 1,  rel_heights = c(.1,.9))
p12
save_plot(p12, file="confusion.png", base_width = 14, base_height = 10)

truth_table_p_0.1

accuracies <- c(95/100, 89/100, 88/100, 80/100, 73/100, 67/100, 60/100, 55/100, 51/100, 44/100)



p_outliers_list <- c(seq(0.1, 1, 0.1))
accuracies <- data.frame(accuracies = accuracies, p_outliers = p_outliers_list)
p1 <- ggplot(data=accuracies) + geom_line(aes(x=p_outliers, y=accuracies), alpha=0.5) +
  labs(x="p_outliers", y="accuracy", title="accuracy across proportion of outliers")
p1
save_plot(p1, file="accuracy.png", base_width = 14, base_height = 10)
#precision = TP/(TP+FP)
precisions <- c(5/5, 9/9, 1, 1, 1, 1, 1, 1, 1, 1)
#recall = TP/(TP+FN)
recalls <- c(5/10, 9/20, 18/30, 20/40, 23/50, 27/33, 30/70, 35/80, 41/90, 44/100)

f1s <- 2 * ((precisions * recalls)/(precisions + recalls))

f1s <- data.frame(f1s = f1s, p_outliers = p_outliers_list)
library(ggplot2)
p1 <- ggplot(data=f1s) + geom_line(aes(x=p_outliers, y=f1s), alpha=0.5) +
  labs(x="p_outliers", y="f1", title="f1 score across proportion of outliers")
p1
save_plot(p1, file="f1.png", base_width = 14, base_height = 10)
