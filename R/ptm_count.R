suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(prob)
  library(stringr)
})


source("./seurat_RNA.R")
source("./seurat_mutations.R")
source("./cell_cycle.R")

## START ##

path.to.levels <- "/<redacted>/SRCdata/EpSRC/EpSRC_RNA_1_outs/filtered_gene_bc_matrices/genome/"
path.to.mut.matrices <- "/<redacted>/working/vartrix_counts_CTRL/EPSRC1rna/"
suob.rna <- ReadRNA(path.to.levels)
suob.mut <- ReadMut(path.to.mut.matrices)
suob.mut.meta <- read.table("/<redacted>/working/PTM_overlap/EPSRC1_counted.csv", header=TRUE, sep=",")
suob.mut.meta <- indexMutMetaUnder(suob.mut.meta)
suob.mut[["ALT"]] <- suob.mut[["ALT"]][which(rownames(suob.mut[["ALT"]]) %in% suob.mut.meta$index),]
suob.mut[["REF"]] <- suob.mut[["REF"]][which(rownames(suob.mut[["REF"]]) %in% suob.mut.meta$index),]
suob <- loadSeuratWithMutations(suob.rna, suob.mut, "EpSRCrna1")

## RNA PROCESSING ##

suob <- runSeuratRNA(suob, standard_filters=TRUE)

suob <- cellCycleAnnotate(suob)
suob <- regressCellCycle(suob)

## MUTATION PROCESSING

#suob.mut.matched.matched <- suob.mut.matched[,(colnames(suob.mut.matched) %in% colnames(suob))]

#suob[["MUT"]] <- CreateAssayObject(counts = suob.mut.matched.matched)


suob.mut.out <- as.data.frame(suob[["ALT"]]@counts)
suob.mut.final <- as.data.frame(suob.mut.out > 0)
suob.mut.final$Count <- rowSums(suob.mut.final)
suob.mut.final$Mean <- rowMeans(suob.mut.out)
suob.mut.filtered <- suob.mut.final[suob.mut.final$Count >= 10,]
suob.mut.filtered <- subset(suob.mut.filtered, select = -c(Count, Mean))
suob.mut.filtered.t <- as.data.frame(t(suob.mut.filtered))
suob.mut.filtered.t$Clusters <- Idents(suob)

cluster.counts <- aggregate(. ~ Clusters, suob.mut.filtered.t, sum)
cluster.counts <- as.data.frame(cluster.counts)
cluster.counts$Total <- (suob.mut.filtered.t %>% dplyr::filter(!is.na(Clusters)) %>% count(Clusters))$n
p.values <- list()
for (mut in colnames(cluster.counts[-which(names(cluster.counts) %in% c("Total","Columns", "Clusters"))]))
{
  
  p.values[mut] <- chisq.test(data.frame(cluster.counts[which(names(cluster.counts) == mut)],cluster.counts$Total-cluster.counts[which(names(cluster.counts) == mut)]))$p.value
}
p.values <- t(as.data.frame(t(p.values)))
#p.values.filtered <- t(as.data.frame(t(p.values[p.values[,1] <= 10^-10,])))
p.values.filtered <- p.values
colnames(x=p.values.filtered) <- c("alt.p")

suob.mut.ref.out <- subset(t(as.data.frame(suob[["REF"]]@counts)), select=rownames(p.values))
suob.mut.ref.filtered.t <- as.data.frame(suob.mut.ref.out > 0)
suob.mut.ref.filtered.t$Clusters <- Idents(suob)

cluster.counts.ref <- aggregate(. ~ Clusters, suob.mut.ref.filtered.t, sum)
cluster.counts.ref <- as.data.frame(cluster.counts.ref)
p.values.ref <- list()
for (mut in colnames(cluster.counts[-which(names(cluster.counts) %in% c("Total","Columns", "Clusters"))]))
{
  #p.values.ref[mut] <- chisq.test(data.frame(cluster.counts[which(names(cluster.counts) == mut)],cluster.counts.ref[which(names(cluster.counts.ref) == mut)]))$p.value
  p.values.ref[mut] <- fisher.test(cbind(cluster.counts[which(names(cluster.counts) == mut)],cluster.counts.ref[which(names(cluster.counts.ref) == mut)]), workspace=1e+7, hybrid=TRUE)$p.value
}
p.values.ref <- t(as.data.frame(t(p.values.ref)))
colnames(x=p.values.ref) <- c("ref.p")
p.values.ref.filtered <- t(as.data.frame(t(p.values.ref[p.values.ref[,1] <= 0.1,])))
colnames(x=p.values.filtered) <- c("alt.p")

rownames(p.values.ref) <- gsub("-", "_", rownames(p.values.ref))
p.values.ref <- data.frame(p.values.ref)
p.values.ref$gene <- suob.mut.meta.copy[rownames(p.values.ref),]$gene

print(p.values.ref)

suob.mut.ref.out <- subset(t(as.data.frame(suob[["REF"]]@counts)), select=rownames(p.values))
suob.mut.ref.reads <- as.data.frame(suob.mut.ref.out)
suob.mut.ref.reads$Clusters <- Idents(suob)
cluster.counts.ref.reads <- aggregate(. ~ Clusters, suob.mut.ref.reads, sum)
cluster.counts.ref.reads <- as.data.frame(cluster.counts.ref.reads)

suob.mut.out <- subset(t(as.data.frame(suob[["ALT"]]@counts)), select=rownames(p.values))
suob.mut.reads <- as.data.frame(suob.mut.out)
suob.mut.reads$Clusters <- Idents(suob)
cluster.counts.reads <- aggregate(. ~ Clusters, suob.mut.reads, sum)
cluster.counts.reads <- as.data.frame(cluster.counts.reads)

cluster.counts.reads/cluster.counts[-which(names(cluster.counts) %in% c("Total"))]
cluster.counts.ref.reads/cluster.counts.ref

rownames(p.values.ref) <- gsub("-", "_", rownames(p.values.ref))
suob.mut.meta$ref.p <- data.frame(p.values.ref)[rownames(suob.mut[["ALT"]]),"ref.p"]

# TODO fix for case where length is zero, fix for case where length is one
i <- 0
while (i*9 < nrow(p.values.ref.filtered))
{
  j <- min(c((i+1)*9, nrow(p.values.ref.filtered)))
  print(i*9+1)
  print(j)
  X11()
  p = FeaturePlot(suob, features = paste("alt_", rownames(p.values.ref.filtered)[(i*9+1):(j)], sep=""), ncol=3)
  print(p)
  i <- i + 1
}

i <- 0
while (i*9 < nrow(p.values.ref.filtered))
{
  j <- min(c((i+1)*9, nrow(p.values.ref.filtered)))
  print(i*9+1)
  print(j)
  X11()
  p = FeaturePlot(suob, features = paste("ref_", rownames(p.values.ref.filtered)[(i*9+1):(j)], sep=""), ncol=3)
  print(p)
  i <- i + 1
}

# hold up on that filtering, you can also have the opposite direction (ie. low p value on ref, high p-value on alt) so you should be checking on this as well (this requires additional processing to get a separate list of filtered variants)
# consider moving to some form of differential expression, this will allow you to see which clusters are affected

typeline("Done, press enter to close")