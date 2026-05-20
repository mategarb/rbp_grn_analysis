################################################################################
# Consensus RBP regulatory-network analysis
#
# Purpose:
#   Build and compare HepG2 and K562 consensus regulatory networks, visualize
#   network statistics, annotate selected interactions with eCLIP/RAP-seq data,
#   and perform downstream validation in LIHC/GTEx datasets.
#
# Notes for reuse:
#   1. Update only the paths and parameters in the Configuration section.
#   2. All user-specific absolute paths were replaced with dummy paths.
#   3. The analysis logic from the original script was preserved as much as
#      possible; commented-out exploratory lines were left in place when useful.
################################################################################

# Required packages ------------------------------------------------------------
required_packages <- c(
  "R.matlab",
  "pheatmap",
  "tidyverse",
  "igraph",
  "reshape2",
  "ape",
  "venn",
  "ggvenn",
  "arules",
  "pracma",
  "e1071",
  "VennDiagram",
  "gprofiler2",
  "ggplot2",
  "ggplotify",
  "ggraph",
  "Ckmeans.1d.dp",
  "data.table",
  "rtracklayer",
  "AnnotationHub",
  "VariantAnnotation",
  "GenomicRanges",
  "RCAS",
  "readxl",
  "ggpubr",
  "eulerr",
  "ggsci",
  "ggnet",
  "network",
  "sna",
  "ggnetwork",
  "IRanges",
  "Gviz",
  "trackViewer",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "svglite",
  "survival",
  "survminer",
  "scales",
  "classInt",
  "caTools",
  "party",
  "dplyr",
  "magrittr",
  "ISLR",
  "rpart",
  "rpart.plot",
  "tidyr",
  "gridExtra",
  "grid",
  "lattice",
  "caret",
  "mclust",
  "DOSE",
  "enrichplot",
  "AnnotationDbi",
  "harmonicmeanp",
  "ggpmisc",
  "plotly",
  "ggrepel",
  "clusterProfiler",
  "stringr",
  "arcdiagram",
  "rrvgo"
)

invisible(lapply(required_packages, require, character.only = TRUE))

# Configuration ----------------------------------------------------------------
# Replace the dummy paths below with paths on your machine before running.
paths <- list(
  project_dir = "/path/to/2p_ENCODE_project",
  functions_dir = "/path/to/2p_ENCODE_project/scripts/funs",
  output_dir = "/path/to/output",
  
  # Core input files
  tf_names_file = "/path/to/2p_ENCODE_project/main_files_results_ENCODE/some_files/TF_names_v_1.01.txt",
  hepg2_gene_names_mat = "/path/to/2p_ENCODE_project/main_files_results_ENCODE/networks_encode3/geneNames_hepg2_prank.mat",
  k562_gene_names_mat = "/path/to/2p_ENCODE_project/main_files_results_ENCODE/networks_encode3/geneNames_k562_prank.mat",
  hepg2_expression_file = "/path/to/2p_ENCODE_project/datasets/ymatrix_hepg2.csv",
  k562_expression_file = "/path/to/2p_ENCODE_project/datasets/ymatrix_k562.csv",
  hepg2_network_dir = "/path/to/2p_ENCODE_project/ENCODE_realNets/HepG2",
  k562_network_dir = "/path/to/2p_ENCODE_project/ENCODE_realNets/K562",
  
  # eCLIP/RAP-seq and annotation files
  eclip_dir = "/path/to/2p_ENCODE_project/datasets/eCLIP_data",
  eclip_metadata_file = "/path/to/2p_ENCODE_project/datasets/eCLIP_data/metadata.tsv",
  rapseq_peaks_dir = "/path/to/RAP-seq_riccardo/peaks",
  gencode_gtf_file = "/path/to/2p_ENCODE_project/datasets/gencode.v45.annotation.gtf.gz",
  rbp_list_file = "/path/to/2p_ENCODE_project/datasets/rbps_info/210329_Table_S1_hRBP_list.xlsx",
  
  # Validation input files
  gtex_pheno_file = "/path/to/2p_ENCODE_project/datasets/validation_UCSC_LIHC_GTEx/GTEX_phenotype.gz",
  lihc_survival_file = "/path/to/2p_ENCODE_project/datasets/validation_UCSC_LIHC_GTEx/survival_LIHC_survival.txt",
  gene_probe_map_file = "/path/to/2p_ENCODE_project/datasets/validation_UCSC_LIHC_GTEx/probeMap_gencode.v23.annotation.gene.probemap",
  liver_tcga_gtex_rds = "/path/to/2p_ENCODE_project/results/liver_tcga_gtex.rds",
  all_interactions_rds = "/path/to/2p_ENCODE_project/results/all_inters_hepg2.rds",
  cancer_literature_rbp_file = "/path/to/2p_ENCODE_project/datasets/rbps_info/TheNumberOfCancerRelevantLiteraturesOfAllRBPs.csv",
  cancer_de_rbp_file = "/path/to/2p_ENCODE_project/datasets/rbps_info/TheNumberOfDifferentiallyExpressedRBPsOfAllCancerTypes.csv",
  functional_network_file = "/path/to/2p_ENCODE_project/datasets/FC5.0_H.sapiens_full.gz",
  msigdb_gmt_file = "/path/to/2p_ENCODE_project/MSigDB/c2.cp.v2023.2.Hs.symbols.gmt"
)

params <- list(
  consensus_network_to_plot = 3,
  hepg2_network_threshold = 0.3,
  number_of_expression_replicates = 232,
  fdr_cutoff = 0.05
)

dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$output_dir, "supplementary_tables"), recursive = TRUE, showWarnings = FALSE)

## edge with eCLIP
# Helper functions ------------------------------------------------------------
# convert networks into directed and discretized
catnet <- function(net){
  net[net < 0] <- -1
  net[net > 0] <- 1
  return(net)}

source(file.path(paths$functions_dir, "intersectNets.R"))
source(file.path(paths$functions_dir, "intersectVarSizeNets.R"))
source(file.path(paths$functions_dir, "hypergeomVenn.R"))
# Load external helper scripts ------------------------------------------------

# read list of transcription factors
TFs <- read.table(paths$tf_names_file)
TFs <- as.character(as.matrix(TFs)) #http://humantfs.ccbr.utoronto.ca/

# read gene names
geneNames_hepg2 <- readMat(paths$hepg2_gene_names_mat)
geneNames_hepg2 <- trimws(geneNames_hepg2$nams.hepg2) #%>% as.vector()
geneNames_k562 <- readMat(paths$k562_gene_names_mat)
geneNames_k562 <- trimws(geneNames_k562$nams.k562)


# read HepG2 data
hepg2 <- read.csv2(file=paths$hepg2_expression_file, header=T, sep="\t",row.names = 1)
hepg2 <- as.data.frame(sapply(hepg2, as.numeric))
#rownames(hepg2) <- rownames(hepg20)
#gensToRem1 <- setdiff(rownames(hepg2),geneNames_hepg2)
#rtorem <- match(gensToRem1,rownames(hepg2))
#ctorem <- c(match(gensToRem1,rownames(hepg2)), match(gensToRem1,rownames(hepg2))+length(rownames(hepg2)))
#hepg2 <- hepg2[-rtorem,-ctorem]
#Ph <- -cbind(eye(dim(hepg2)[1]),eye(dim(hepg2)[1]))
#Ah <- -Ph %*% pinv(hepg2 %>% as.matrix) # infering sign for networks based on ls
# read K562 data
k562 <- read.csv2(file=paths$k562_expression_file, header=T, sep="\t",row.names = 1)
k562 <- as.data.frame(sapply(k562, as.numeric))  
#rownames(k562) <- rownames(k5620)
#gensToRem2 <- setdiff(rownames(k562),geneNames_k562)
#rtorem <- match(gensToRem2,rownames(k562))
#ctorem <- c(match(gensToRem2,rownames(k562)), match(gensToRem2,rownames(k562))+length(rownames(k562)))
#k562 <- k562[-rtorem,-ctorem] 
#Pk <- -cbind(eye(dim(k562)[1]),eye(dim(k562)[1])) # infering sign for networks based on ls
#Ak <- -Pk %*% pinv(k562 %>% as.matrix)


# common part of genes between cell lines
x <- list(HepG2=geneNames_hepg2, K562=geneNames_k562)
p01 <- ggvenn(x, fill_color = c("#e39384", "#6fc5c9"), fill_alpha = 0.5,
              stroke_size = 0.5, set_name_size = 4) + theme(text = element_text(size = 16))



############################
# create inputs from HepG2 #
############################
### read HepG2 networks

path <- paths$hepg2_network_dir
allfs <- list.files(path)
namns <- lapply(strsplit(allfs,"_"), function(x) x[2]) %>% unlist

lnets <- list()
lnet <- list()
for(i in 1:length(namns))
{
  ptmp <- paste0(path, allfs[i])
  
  tmpnet <- read.table(ptmp, sep=",")
  lnet[[i]] <- catnet(tmpnet)
  lnets[[i]] <- tmpnet
}


dfh10 <- lapply(lnet, function(x) as.numeric(as.logical(as.vector(as.matrix(x)))))
dfh10 <- do.call(data.frame, dfh10)
colnames(dfh10) <- namns

dfl10 <- lapply(lnet, function(x) as.vector(as.matrix(x)))
dfl10 <- do.call(data.frame, dfl10)
colnames(dfl10) <- namns

dfll10 <- lnet
names(dfll10) <- namns

###########################
# create inputs from K562 #
###########################

pathk <- paths$k562_network_dir
allfsk <- list.files(pathk)
namnsk <- lapply(strsplit(allfsk,"_"), function(x) x[2]) %>% unlist

lnetk <- list()
for(i in 1:length(namnsk))
{
  ptmpk <- paste0(pathk, allfsk[i])
  
  tmpnetk <- read.table(ptmpk, sep=",")
  lnetk[[i]] <- catnet(tmpnetk)
}


dfh10k <- lapply(lnetk, function(x) as.numeric(as.logical(as.vector(as.matrix(x)))))
dfh10k <- do.call(data.frame, dfh10k)
colnames(dfh10k) <- namnsk

dfl10k <- lapply(lnetk, function(x) as.vector(as.matrix(x)))
dfl10k <- do.call(data.frame, dfl10k)
colnames(dfl10k) <- namnsk

dfll10k <- lnetk
names(dfll10k) <- namnsk
############################

# hierarchical clustering in order to find groups of similar methods
distance_math <- hamming.distance(dfh10 %>% t) %>% as.dist
Hierar_cl <- hclust(distance_math, method = "ward.D")

hadata <- hamming.distance(dfh10 %>% t)
dstc <- hadata[upper.tri(hadata)]
kmod <- Ckmeans.1d.dp(dstc, c(3,9))

c1 <- dstc[kmod$cluster %in% which(kmod$centers<mean(dstc)-2*sd(dstc))]


labelso <- hamming.distance(dfh10 %>% t)
labelso[which(hadata %in% c1)] <- "•"
labelso[which(labelso != "•")] <- ""
p11 <- pheatmap(hadata,  clustering_method="ward.D",fontsize = 14,
                color=colorRampPalette(c("#a561b0", "white", "#58ba54"))(20), border_color = "black",
                number_color = "#FF8000", fontsize_number = 50) %>% as.ggplot()
p11 <- p11 + theme(text = element_text(size = 16))
##k562
distance_math <- hamming.distance(dfh10k %>% t) %>% as.dist
Hierar_cl <- hclust(distance_math, method = "ward.D")
kadata <- hamming.distance(dfh10k %>% t)
p22 <- pheatmap(kadata,  clustering_method="ward.D",fontsize = 14,
                color=colorRampPalette(c("#a561b0", "white", "#58ba54"))(20), border_color = "black",
                number_color = "#FF8000", fontsize_number = 50) %>% as.ggplot()
p22 <- p22 + theme(text = element_text(size = 16))


# HepG2
netsi2 <- lapply(dfll10, abs)
avgNets <- Reduce('+', netsi2)/length(netsi2) # avg. occurrence of a link (frequency of a link)
allthrs <- as.numeric(names(table(sort(unique(abs(avgNets[avgNets!=0])))))) # kind of zetas
inets_hepg2 <- list()
for(i in 1:length(allthrs)){
  
  tmpAvgNets <- avgNets
  tmpAvgNets[abs(tmpAvgNets) < allthrs[i]] <- 0 # keep at least 11, keep at least 10, etc.
  
  tmpAvgNets[tmpAvgNets < 0] <- 1
  tmpAvgNets[tmpAvgNets > 0] <- 1
  
  inets_hepg2[[i]] <- tmpAvgNets
}

S_all <- as.numeric(lapply(inets_hepg2,function(x) mean(apply(x,2,function(x) sum(x!=0))))) %>% round(digits = 2)

lapply(lapply(lapply(inets_hepg2, function(x) apply(x,2,function(x) sum(x!=0))), function(x) x[x!=0]),mean)

lapply(lapply(inets_hepg2, function(x) apply(x,2,function(x) sum(x!=0))), sum)
netc <- inets_hepg2[[5]]
rownames(netc) <- geneNames_hepg2
colnames(netc) <- geneNames_hepg2
g <- graph.adjacency(netc %>% t %>% as.matrix)
alin <- get.edgelist(g)

table(lapply(dfll10, function(x) x[which(geneNames_hepg2=="PES1"), which(geneNames_hepg2=="AQR")]) %>% unlist)
table(lapply(dfll10, function(x) x[which(geneNames_hepg2=="KIF1C"), which(geneNames_hepg2=="RBM39")]) %>% unlist)
table(lapply(dfll10, function(x) x[which(geneNames_hepg2=="FASTKD1"), which(geneNames_hepg2=="RPS10")]) %>% unlist)


# sign for the cGRN
hnet <- avgNets %>% as.matrix #%>% t
#hnet[hnet < 0.1] <- 0 
rownames(hnet) <- geneNames_hepg2
colnames(hnet) <- geneNames_hepg2

hexprs <- read.table(paths$hepg2_expression_file)

hexprs1 <- (hexprs[,1:232] + hexprs[,233:464])/2


## sign of a link
alledgs <- which(hnet!=0, arr.ind = T)
#inde <- hnet[alledgs[1,1], alledgs[1,2]]
#nmne <- c(geneNames_hepg2[alledgs[1,1]], geneNames_hepg2[alledgs[1,2]])

segds <- list()
tmpsign <- c()
pcr <- c()
scr <- c()
for(i in 1:dim(alledgs)[1]){
  for(j in 1:length(dfll10)){
    tmpnet <- dfll10[[j]]
    rownames(tmpnet) <- geneNames_hepg2
    colnames(tmpnet) <- geneNames_hepg2
    tmpsign[j] <- tmpnet[alledgs[i,1], alledgs[i,2]]
  }
  segds[[i]] <- tmpsign
  
  tmpgA <- hexprs1[alledgs[i,1],] %>% as.matrix %>% as.numeric
  tmpgB <- hexprs1[alledgs[i,2],] %>% as.matrix %>% as.numeric
  pcr[i] <- cor(tmpgA, tmpgB, method = 'pearson')
  scr[i] <- cor(tmpgA, tmpgB, method = 'spearman')
}
msedgs <- do.call(rbind,segds)
msedgs2 <- msedgs[,-which(namns %in% c("CART","neunetreg"))]
msedgs2 <- cbind(msedgs2, scr %>% catnet %>% as.data.frame)

tmpscr <- scr %>% catnet
sgn <- c()
for(i in 1:dim(alledgs)[1]){
  tmpt <- table(msedgs2[i,] %>% as.matrix)
  tmpt <- tmpt[-which(names(tmpt) == 0)]
  
  if(length(names(tmpt)) > 1){
    if(tmpt[which(names(tmpt)=="-1")] == tmpt[which(names(tmpt)=="1")]){
      sgn[i] <- tmpscr[i]
    }else if(tmpt[which(names(tmpt)=="-1")] < tmpt[which(names(tmpt)=="1")]){
      sgn[i] <- 1
    }else{
      sgn[i] <- -1
    }
  }else{
    sgn[i] <- names(tmpt) %>% as.numeric
  }
}

for(i in 1:dim(alledgs)[1]){
  hnet[alledgs[i,1], alledgs[i,2]] <- hnet[alledgs[i,1], alledgs[i,2]]*sgn[i]
}

#write.csv2(hnet, file = file.path(paths$output_dir, "Tab.S1_adjacencymatrix_hepg2_signed.csv"))

hnet2 <- hnet %>% t
hnet2[abs(hnet2) < 0.5] <- 0
#write.csv(hnet2, file = file.path(paths$output_dir, "Tab.S1_adjacencymatrix_hepg2_signed_t_cluereg.csv"), quote = F)


## K562

netsi2k <- lapply(dfll10k, abs)
avgNetsk <- Reduce('+', netsi2k)/length(netsi2k) # avg. occurrence of a link (frequency of a link)
allthrsk <- as.numeric(names(table(sort(unique(abs(avgNetsk[avgNetsk!=0])))))) # kind of zetas
inets_k562 <- list()
for(i in 1:length(allthrs)){
  
  tmpAvgNets <- avgNetsk
  tmpAvgNets[abs(tmpAvgNets) < allthrsk[i]] <- 0 # keep at least 11, keep at least 10, etc.
  
  tmpAvgNets[tmpAvgNets < 0] <- 1
  tmpAvgNets[tmpAvgNets > 0] <- 1
  
  inets_k562[[i]] <- tmpAvgNets
}

S_allk <- as.numeric(lapply(inets_k562,function(x) mean(apply(x,2,function(x) sum(x!=0))))) %>% round(digits = 2)
netc <- inets_k562[[7]]
rownames(netc) <- geneNames_k562 
colnames(netc) <- geneNames_k562 
g <- graph.adjacency(netc %>% as.matrix)
get.edgelist(g)

# sign
#netsign <- Reduce('+', lnets)
#rownames(netsign) <- geneNames_hepg2
#colnames(netsign) <- geneNames_hepg2


# save networks to files
rownames(avgNets) <- geneNames_hepg2
colnames(avgNets) <- geneNames_hepg2

#write.csv2(avgNets, file = file.path(paths$output_dir, "supplementary_tables/Tab.S1_adjacencymatrix_hepg2.csv"))

rownames(avgNetsk) <- geneNames_k562 
colnames(avgNetsk) <- geneNames_k562 

#write.csv2(avgNetsk, file = file.path(paths$output_dir, "supplementary_tables/Tab.S2_adjacencymatrix_k562.csv"))

# sign for the cGRN - k562
knet <- avgNetsk %>% as.matrix #%>% t
#hnet[hnet < 0.1] <- 0 
rownames(knet) <- geneNames_k562
colnames(knet) <- geneNames_k562

kexprs <- read.table(paths$k562_expression_file)

kexprs1 <- (kexprs[,1:232] + kexprs[,233:464])/2


## sign of a link
alledgsk <- which(knet!=0, arr.ind = T)
#inde <- hnet[alledgs[1,1], alledgs[1,2]]
#nmne <- c(geneNames_hepg2[alledgs[1,1]], geneNames_hepg2[alledgs[1,2]])

segds <- list()
tmpsign <- c()
pcr <- c()
scr <- c()
for(i in 1:dim(alledgsk)[1]){
  for(j in 1:length(dfll10k)){
    tmpnet <- dfll10k[[j]]
    rownames(tmpnet) <- geneNames_k562
    colnames(tmpnet) <- geneNames_k562
    tmpsign[j] <- tmpnet[alledgsk[i,1], alledgsk[i,2]]
  }
  segds[[i]] <- tmpsign
  
  tmpgA <- kexprs1[alledgsk[i,1],] %>% as.matrix %>% as.numeric
  tmpgB <- kexprs1[alledgsk[i,2],] %>% as.matrix %>% as.numeric
  pcr[i] <- cor(tmpgA, tmpgB, method = 'pearson')
  scr[i] <- cor(tmpgA, tmpgB, method = 'spearman')
}
msedgs <- do.call(rbind,segds)
msedgs2 <- msedgs[,-which(namns %in% c("CART","neunetreg"))]
msedgs2 <- cbind(msedgs2, scr %>% catnet %>% as.data.frame)

tmpscr <- scr %>% catnet
sgn <- c()
for(i in 1:dim(alledgsk)[1]){
  tmpt <- table(msedgs2[i,] %>% as.matrix)
  tmpt <- tmpt[-which(names(tmpt) == 0)]
  
  if(length(names(tmpt)) > 1){
    if(tmpt[which(names(tmpt)=="-1")] == tmpt[which(names(tmpt)=="1")]){
      sgn[i] <- tmpscr[i]
    }else if(tmpt[which(names(tmpt)=="-1")] < tmpt[which(names(tmpt)=="1")]){
      sgn[i] <- 1
    }else{
      sgn[i] <- -1
    }
  }else{
    sgn[i] <- names(tmpt) %>% as.numeric
  }
}

for(i in 1:dim(alledgsk)[1]){
  knet[alledgsk[i,1], alledgsk[i,2]] <- knet[alledgsk[i,1], alledgsk[i,2]]*sgn[i]
}
write.csv2(knet, file = file.path(paths$output_dir, "Tab.S1_adjacencymatrix_k562_signed.csv"))

###########################################
#### intersect HepG2 and K562 networks ####
###########################################
wnet <- params$consensus_network_to_plot # links must be supported by at least this many methods
ih <- inets_hepg2[[wnet]]
ik <- inets_k562[[wnet]]
colnames(ih) <- geneNames_hepg2
rownames(ih) <- geneNames_hepg2
colnames(ik) <- geneNames_k562
rownames(ik) <- geneNames_k562
inter_hk <- intersectVarSizeNets(ih %>% t, ik %>% t)

hepg2 <- lapply(inets_hepg2, function(x) mean(colSums(catnet(x)))) %>% as.numeric %>% round(digits = 5)
k562 <- lapply(inets_k562, function(x) mean(colSums(catnet(x)))) %>% as.numeric %>% round(digits = 5)
df2 <- data.frame(HepG2 = hepg2, K562 = k562, nets=as.character(1:length(hepg2))) %>% melt(id="nets")
p2 <- ggplot(df2, aes(x=nets, y=value, fill = variable))+
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#e39384", "#6fc5c9"), name="") +
  geom_text(aes(label = round(value, digits=2)),position = position_dodge(width = 0.9), vjust = -0.25)+
  ylab("node degree")+
  xlab("number of consensus links")+ theme_pubclean()+
  theme(legend.title=element_blank(), text = element_text(size = 16), legend.position="top")
a1 <- unname(table(catnet(inets_hepg2[[wnet]]) %>% as.matrix)[2])
a2 <- unname(table(catnet(inets_k562[[wnet]]) %>% as.matrix)[2])
ca1 <- unname(table(inter_hk %>% as.matrix)[names(table(inter_hk %>% as.matrix)) == "common_link"])
ca2 <- unname(table(inter_hk %>% as.matrix)[names(table(inter_hk %>% as.matrix)) == "common_link"]) + unname(table(inter_hk %>% as.matrix)[names(table(inter_hk %>% as.matrix)) == "uncommon_link"])

vd <- euler(c(HepG2 = a1,          # Draw pairwise venn diagram
              K562 = a2,
              "HepG2&K562" = ca1))

p3 <- plot(vd, edges=c("black","black"), fills=list(fill=c("#e39384", "#6fc5c9"), alpha=0.5),
           quantities = list(type = c("counts","percent")))

vd2 <- euler(c(HepG2 = a1,          # Draw pairwise venn diagram
               K562 = a2,
               "HepG2&K562" = ca2))

p4 <- plot(vd2, fills=c("#674C47", "#880808"),quantities = list(type = c("counts","percent")), alpha = 0.5)


# common and uncommon links content
val3 <- table(inter_hk, useNA = "always")
val4 <- melt(val3)

p5 <- ggplot(val4[c(-1,-5),], aes(x="", y=value, fill=inter_hk)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+ 
  theme_void() +
  scale_fill_manual(values=c("#8c200f","#7a6663","#d9553f"),name = "link", labels = c("common", "unknown","uncommon")) + 
  geom_text(aes(label = paste(round(val4[c(-1,-5),]$value/sum(val4[c(-1,-5),]$value)*100, digits=1), "%"), colour="white"), position = position_stack(vjust = 0.5), colour="white", size=5)

## how methods contributed to a given link inference, only TPs
permetl <- list()
for(j in 1:length(inets_hepg2)){
  permet <- c()
  for(i in 1:length(lnet)){
    permet[i] <- length(which((inets_hepg2[[j]] %>% as.matrix %>% as.vector %>% abs)==1 & (lnet[[i]] %>% as.matrix %>% as.vector %>% abs)==1))
  }
  names(permet) <- namns
  permetl[[j]] <- permet/sum(permet)
}

df <- do.call(rbind,permetl) %>% as.data.frame
df$id <- paste0(1:8,"/10")
df2 <- melt(df)


#df$subject <- factor(df$subject) 
#df$credit <- factor(df$credit)  
cls <- c("#d65159", "#1dd582", "#628c9c", "#e093b4", "#b7b62f", "#051453", "#f8ffa3","#D55E00", "#94c6ff", "#D318CA", "#5d8551","#5c2e2e" )

ppie <- ggplot(data=df2, aes(x=" ", y=value, group=variable, colour=variable, fill=variable)) + 
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +  
  facet_wrap(.~ id, ncol = 3) + 
  scale_fill_manual(values = cls, name="method") +
  scale_colour_manual(values = rep("black",length(cls)), name="method")+theme_void() +
  theme(text = element_text(size=16))

# k562
permetlk <- list()
for(j in 1:length(inets_k562)){
  permetk <- c()
  for(i in 1:length(lnet)){
    permetk[i] <- length(which((inets_k562[[j]] %>% as.matrix %>% as.vector %>% abs)==1 & (lnetk[[i]] %>% as.matrix %>% as.vector %>% abs)==1))
  }
  names(permetk) <- namns
  permetlk[[j]] <- permetk/sum(permetk)
}

df <- do.call(rbind,permetlk) %>% as.data.frame
df$id <- paste0(1:8,"/10")
df2 <- melt(df)

ppie2 <- ggplot(data=df2, aes(x=" ", y=value, group=variable, colour=variable, fill=variable)) + 
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +  
  facet_wrap(.~ id, ncol = 3) + 
  scale_fill_manual(values = cls, name="method") +
  scale_colour_manual(values = rep("black",length(cls)), name="method")+theme_void() +
  theme(text = element_text(size=16))
###########

gup <- ggarrange(p11, ppie,
                 labels = c("A", "B"),
                 ncol = 2, nrow = 1)

gmid <- ggarrange(p22, ppie2,
                  labels = c("C", "D"),
                  ncol = 2, nrow = 1)

gdo <- ggarrange(ggarrange(p01, p3,
                           labels = c("E","F"),
                           ncol = 1, nrow = 2), p2,
                 labels = c("","G"),
                 ncol = 2, nrow = 1, widths = c(0.7,1.3))

ggall <- ggarrange(gup,gmid, gdo, nrow=3,ncol=1)
ggsave(ggall, file = file.path(paths$output_dir, "figures/encode_stats_realNets.svg"), limitsize = T, width = 12, height = 14, dpi =700, bg = 'white')

##################################################
######### NETWORK VISUALIZATION ##################
##################################################

# give a sign to the HepG2 network
#signsNet <- catnet(Ah)
#displayNet <- abs(inets_hepg2[[2]]) %>% as.matrix
#displayNet <- displayNet * signsNet

#newVals_uncommon <- displayNet[inter_hk == "uncommon_link"]
#newVals_uncommon[newVals_uncommon < 0] <- newVals_uncommon[newVals_uncommon < 0]-5
#newVals_uncommon[newVals_uncommon > 0] <- newVals_uncommon[newVals_uncommon > 0]+5

#newVals_common <- displayNet[inter_hk == "common_link"]
#newVals_common[newVals_common < 0] <- newVals_common[newVals_common < 0]-10
#newVals_common[newVals_common > 0] <- newVals_common[newVals_common > 0]+10

#displayNet0 <- zeros(size(displayNet)[1])
#displayNet0[inter_hk == "uncommon_link"] <- newVals_uncommon
#displayNet0[inter_hk == "common_link"] <- newVals_common
#displayNet0[inter_hk == "uncommon_link"] <- 1
#displayNet0[inter_hk == "common_link"] <- 2
#colnames(displayNet0) <- colnames(displayNet)
#rownames(displayNet0) <- rownames(displayNet)
# display the network of link frequencies
#pheatmap(displayNet, cluster_rows=FALSE, cluster_cols=FALSE)

#################################################
####### displaying a network in cytoscape #######
#################################################
hnet <- avgNets %>% as.matrix %>% t
hnet[hnet < params$hepg2_network_threshold] <- 0 
rownames(hnet) <- geneNames_hepg2
colnames(hnet) <- geneNames_hepg2
hig <- hnet %>% 
  graph_from_adjacency_matrix(mode = "directed", weighted = TRUE)

Isolated = which(igraph::degree(hig)==0)
hig = igraph::delete.vertices(hig, Isolated)

alledgs <- which(hnet!=0, arr.ind = T)

# remove genes without any links
nrem <- intersect(unname(which(colMeans(hnet)==0)), unname(which(rowMeans(hnet)==0)))
hnet <- hnet[-nrem,-nrem]
hexprs <- read.table(paths$hepg2_expression_file)

hexprs1 <- (hexprs[,1:232] + hexprs[,233:464])/2

V(hig)$Node_degree = igraph::degree(hig)
mexps <- apply(hexprs, 1, median, na.rm=T)
mexps <- mexps[match(names(V(hig)), names(mexps))] %>% unname
V(hig)$mfc = mexps
#colfunc <- colorRampPalette(c("#C19DEE", "#FBAA11"))

vxnms2 <- names(V(hig))
tfrbp <- rep("RBP",length(vxnms2))
tfrbp[which(!is.na(match(vxnms2, TFs)))] <- "RBP and TF"
V(hig)$Node_shape = tfrbp

## sign of a link
alledgs <- get.edgelist(hig)
segds <- list()
tmpsign <- c()
pcr <- c()
scr <- c()
for(i in 1:dim(alledgs)[1]){
  for(j in 1:length(dfll10)){
    tmpnet <- dfll10[[j]] %>% t
    rownames(tmpnet) <- geneNames_hepg2
    colnames(tmpnet) <- geneNames_hepg2
    tmpsign[j] <- tmpnet[which(rownames(tmpnet) == alledgs[i,][1]),
                         which(colnames(tmpnet) == alledgs[i,][2])]
  }
  segds[[i]] <- tmpsign
  
  tmpgA <- hexprs1[which(rownames(hexprs1) == alledgs[i,][1]),] %>% as.matrix %>% as.numeric
  tmpgB <- hexprs1[which(rownames(hexprs1) == alledgs[i,][2]),] %>% as.matrix %>% as.numeric
  pcr[i] <- cor(tmpgA, tmpgB, method = 'pearson')
  scr[i] <- cor(tmpgA, tmpgB, method = 'spearman')
}
msedgs <- do.call(rbind,segds)
msedgs2 <- msedgs[,-which(namns %in% c("CART","neunetreg"))]
msedgs2 <- cbind(msedgs2, scr %>% catnet %>% as.data.frame)

tmpscr <- scr %>% catnet
sgn <- c()
for(i in 1:dim(alledgs)[1]){
  tmpt <- table(msedgs2[i,] %>% as.matrix)
  tmpt <- tmpt[-which(names(tmpt) == 0)]
  
  if(length(tmpt)==1){
    sgn[i] <- tmpscr[i]
  }
  else{
    
    if(tmpt[which(names(tmpt)=="-1")] == tmpt[which(names(tmpt)=="1")]){
      sgn[i] <- tmpscr[i]
    }else if(tmpt[which(names(tmpt)=="-1")] < tmpt[which(names(tmpt)=="1")]){
      sgn[i] <- 1
    }else{
      sgn[i] <- -1
    }
  }
}
sgn <- sgn %>% as.character %>% as.factor
levels(sgn) <- c("inhibition","activation")
E(hig)$sign <- sgn


# metadata
meta <- read_tsv(paths$eclip_metadata_file)
#list of RBPs
clsupp <- read_excel(paths$rbp_list_file)
annots <- import.gff(paths$gencode_gtf_file)

exps <- unique(gsub("-human","",meta$`Experiment target`))

ugns <- unique(alledgs[,1])
ugns2 <- ugns[na.omit(match(exps, ugns))]
ogns <- list()

## eclip-seq
for(i in 1:length(ugns2)){
  targ_gen <- ugns2[i]
  meta_sub <- meta[gsub("-human","",meta$`Experiment target`) == targ_gen,]
  
  file <- meta_sub$`File accession`[which(meta_sub$`File format`=="bigBed narrowPeak")]
  countData_bb <- import.bb(paste0(file.path(paths$eclip_dir, ""),file,".bigBed"))
  names(mcols(countData_bb)) <- c("name","score","signalValue","pValue","qValue","peak")
  
  overlaps <- as.data.table(queryGff(queryRegions = countData_bb, gffData = annots))
  overlaps <- overlaps[which(overlaps$gene_type=="protein_coding"),]
  overlaps <- overlaps[which(overlaps$type=="gene"),]
  overlaps <- overlaps[which(as.numeric(p.adjust(10^(-as.numeric(overlaps$query_pValue)), method = "fdr")) <= 0.05),]
  
  targs <- alledgs[which(alledgs[,1] == ugns2[i]),2]
  lbs <- rep("•",length(targs))
  lbs[which(match(targs, overlaps$gene_name)!="NA")] <- "★"
  ogns[[i]] <- lbs
  print(i)
}
names(ogns) <- ugns2

eclipv <- c()
for(i in 1:length(ugns)){
  if(any(ugns2 == ugns[i])){
    eclipv <- c(eclipv, ogns[which(names(ogns) == ugns[i])] %>% unname %>% unlist) 
  }else{
    eclipv <- c(eclipv, c(rep("", length(which(alledgs[,1] == ugns[i])))))
  }
}

## rapseq
raps <- list.files(paths$rapseq_peaks_dir)
allrs2 <- gsub("\\..*","",as.character(as.matrix(raps)))
gnsr <- intersect(allrs2, ugns)
rogns <- list()
for(i in 1:length(gnsr)){
  
  file <- paste0(file.path(paths$rapseq_peaks_dir, ""),gnsr[i], ".peaks.txt")
  overlaps <- read.table(file, header = T)
  overlaps <- overlaps[which(overlaps$gene_type=="protein_coding"),]
  
  targs <- alledgs[which(alledgs[,1] == gnsr[i]),2]
  lbs <- rep("○",length(targs))
  lbs[which(match(targs, overlaps$gene_name)!="NA")] <- "■"
  rogns[[i]] <- lbs
}
names(rogns) <- gnsr

rapv <- c()
for(i in 1:length(ugns)){
  if(any(gnsr == ugns[i])){
    rapv <- c(rapv, rogns[which(names(rogns) == ugns[i])] %>% unname %>% unlist) 
  }else{
    rapv <- c(rapv, c(rep("", length(which(alledgs[,1] == ugns[i])))))
  }
}

E(hig)$eclip <- paste0(eclipv, rapv)
g1 <- ggraph(hig, layout = 'fr', circular=F) +
  geom_edge_arc(aes(width = weight*10, color = sign, label = eclip), label_colour = "#3F3F3F", label_size = 6, 
                alpha = 0.6, strength = 0.3, angle_calc = 'along', label_dodge = unit(2.5, 'mm'),
                arrow = arrow(length = unit(6, 'pt')), end_cap = circle(6, 'pt')) + 
  geom_node_point(aes(size = Node_degree, color = mfc, shape = Node_shape)) + 
  geom_node_text(aes(label = name, size = 2), repel = TRUE, point.padding = unit(0.2, "lines")) +
  theme_graph(base_size=14)+ 
  scale_size_continuous(range = c(3, 6)) +
  scale_colour_gradient2(name="Fold-change",
                         low = "#AA29EA",
                         mid = "#AAAAAA",
                         high = "#FF7E0C")+ 
  scale_edge_width(range = c(0.75, 2.25), name="Cons. links") +
  scale_edge_color_manual(values=c("#FF320D","#0D91FF"), name="Interaction") +
  scale_label_size(range = c(1, 1)) +
  labs(size="Node degree", shape="Node shape")

### PLOT 2

hnet <- inets_hepg2[[3]] %>% as.matrix %>% t
rownames(hnet) <- geneNames_hepg2
colnames(hnet) <- geneNames_hepg2
hig3 <- hnet %>% 
  graph_from_adjacency_matrix(mode = "directed", weighted = TRUE)

hubs <- igraph::degree(hig3, mode="all") %>% sort(decreasing = T)
hubs <- hubs[which(hubs>10)]

tophubs <- sort(hubs[na.omit(match(exps, names(hubs)))], decreasing = T)


g2 <- list()
hfull <- c()
hval <- c()
pval <- c()

alledgs3 <- get.edgelist(hig3)
inters_names <- alledgs3 %>% as.data.frame
colnames(inters_names) <- c("regulator", "target")


outd <- c()
for(i in 1:length(tophubs)){
  targ_gen <- tophubs[i] %>% names
  meta_sub <- meta[gsub("-human","",meta$`Experiment target`) == targ_gen,]
  
  file <- meta_sub$`File accession`[which(meta_sub$`File format`=="bigBed narrowPeak")]
  countData_bb <- import.bb(paste0(file.path(paths$eclip_dir, ""),file,".bigBed"))
  names(mcols(countData_bb)) <- c("name","score","signalValue","pValue","qValue","peak")
  
  overlaps <- as.data.table(queryGff(queryRegions = countData_bb, gffData = annots))
  overlaps <- overlaps[which(overlaps$gene_type=="protein_coding"),]
  overlaps <- overlaps[which(overlaps$type=="gene"),]
  overlaps <- overlaps[which(as.numeric(p.adjust(10^(-as.numeric(overlaps$query_pValue)), method = "fdr")) <= 0.05),]
  
  #hubnet <- inters_names[c(which(inters_names$regulator==targ_gen), which(inters_names$target==targ_gen)),]
  outd[i] <- length(which(inters_names$regulator==targ_gen))/dim(inters_names[c(which(inters_names$regulator==targ_gen), which(inters_names$target==targ_gen)),])[1]
  
  targs <- c(inters_names$target[which(inters_names$regulator==targ_gen)], inters_names$regulator[which(inters_names$target==targ_gen)]) %>% unique # all targets of a hub, undirected
  ogns <- intersect(targs, overlaps$gene_name) # overlap
  
  q <- length(ogns) # number of white balls drawn, common part # 
  k <- intersect(overlaps$gene_name, clsupp$gene_name) %>% length #intersect(overlaps$gene_name, clsupp$gene_name) %>% length #length(targs) # total number of balls drawn EXP2
  
  m <- length(targs) # intersect(overlaps$gene_name, clsupp$gene_name) %>% length # total number of white balls in the urn EXP1
  n <- length(clsupp$gene_name) - m #(clsupp$gene_name %>% length) - m # total number of black balls in the urn # tot
  
  if(q == 0){
    pval[i] <- 1
  }else{
    pval[i] <- 1 - phyper(q,m,n,k)
  }
  
  hfull[i] <- length(targs) # all targets of a hub
  hval[i] <- length(ogns) # overlapping genes
  
  print(i)
}


dfall <- data.frame("gnams"=tophubs %>% names,
                    "notval"=hfull-hval,
                    "val"=hval)
dfall2 <- melt(dfall)
dfall2$gnams <- factor(dfall2$gnams,levels = tophubs %>% names)

veci <- rep(as.numeric(round(hval/hfull,digits=2)*100>=50),2)
veci[rep(as.numeric(round(hval/hfull,digits=2)*100>=75),2)==1] <- 2
veci[veci==0] <- ""
veci[veci=="1"] <- "."
veci[veci=="2"] <- "•"

pval2 <- p.adjust(pval, "fdr")
pval2[pval>0.1] <- "ns"
pval2[pval<=0.1 & pval>=0.05] <- "*"
pval2[pval<=0.05 & pval>=0.01] <- "**"
pval2[pval<=0.01] <- "***"

# Create the barplot
g2 <- ggplot(data=dfall2, aes(x=gnams, y=value, fill=variable)) +
  geom_bar(stat="identity",position = "stack", color = "#2A2A2A")+
  geom_text(aes(y=rep(hfull,2), label=rep(paste0(round(hval/hfull,digits=2)*100,"%"),2)), hjust=1.1, vjust = 0.4,
            color="black", size=3.5)+
  geom_text(aes(y=rep(hfull,2), label=rep(pval2,2)), hjust=-0.3,  vjust = 0.5,
            color="#2A2A2A", size=5)+
  geom_text(aes(y=rep(rep(0,length(hfull)),2), label=rep(paste0("(",round(outd, digits = 2)*100, "%)"),2)), hjust=-0.1,  vjust = 0.3,
            color="white", size=3.5)+
  scale_fill_manual(values=c("#9daab0", "#487DA9"), name="Targets shared\n with eCLIP-seq", aesthetics = "fill", labels=c('No', 'Yes'))+
  theme_classic2()+ rotate_x_text(45)+ylab("Number of targets in eCLIP-seq") + xlab("") +
  theme(text = element_text(size = 16),legend.position="top") + coord_flip()

#### RAP-seq
raps <- list.files(paths$rapseq_peaks_dir)
allrs2 <- gsub("\\..*","",as.character(as.matrix(raps)))
gnsr <- c(allrs2, "IGF2BP2")


hnet <- inets_hepg2[[3]] %>% as.matrix %>% t
rownames(hnet) <- geneNames_hepg2
colnames(hnet) <- geneNames_hepg2
hig3 <- hnet %>% 
  graph_from_adjacency_matrix(mode = "directed", weighted = TRUE)

hubs <- igraph::degree(hig3, mode="all") %>% sort(decreasing = T)
hubs <- hubs[which(hubs>10)]

rtophubs <- sort(hubs[na.omit(match(gnsr, names(hubs)))], decreasing = T)


g22 <- list()
rhfull <- c()
rhval <- c()
rpval <- c()

alledgs3 <- get.edgelist(hig3)
inters_names <- alledgs3 %>% as.data.frame
colnames(inters_names) <- c("regulator", "target")


outd <- c()
for(i in 1:length(rtophubs)){
  
  targ_gen <- names(rtophubs)[i]
  file <- paste0(file.path(paths$rapseq_peaks_dir, ""),targ_gen, ".peaks.txt")
  overlaps <- read.table(file, header = T)
  overlaps <- overlaps[which(overlaps$gene_type=="protein_coding"),]
  outd[i] <- length(which(inters_names$regulator==targ_gen))/dim(inters_names[c(which(inters_names$regulator==targ_gen), which(inters_names$target==targ_gen)),])[1]
  
  targs <- c(inters_names$target[which(inters_names$regulator==targ_gen)], inters_names$regulator[which(inters_names$target==targ_gen)]) %>% unique # all targets of a hub, undirected
  ogns <- intersect(targs, overlaps$gene_name) # overlap
  
  q <- length(ogns) # number of white balls drawn, common part # 
  k <- intersect(overlaps$gene_name, clsupp$gene_name) %>% length #intersect(overlaps$gene_name, clsupp$gene_name) %>% length #length(targs) # total number of balls drawn EXP2
  
  m <- length(targs) # intersect(overlaps$gene_name, clsupp$gene_name) %>% length # total number of white balls in the urn EXP1
  n <- length(clsupp$gene_name) - m #(clsupp$gene_name %>% length) - m # total number of black balls in the urn # tot
  
  if(q == 0){
    rpval[i] <- 1
  }else{
    rpval[i] <- 1 - phyper(q,m,n,k)
  }
  
  rhfull[i] <- length(targs)
  rhval[i] <- length(ogns)
  
  print(i)
}


rdfall <- data.frame("gnams"=rtophubs %>% names,
                     "notval"=rhfull-rhval,
                     "val"=rhval)
rdfall2 <- melt(rdfall)
rdfall2$gnams <- factor(rdfall2$gnams,levels = rtophubs %>% names)

veci <- rep(as.numeric(round(rhval/rhfull,digits=2)*100>=50),2)
veci[rep(as.numeric(round(rhval/rhfull,digits=2)*100>=75),2)==1] <- 2
veci[veci==0] <- ""
veci[veci=="1"] <- "."
veci[veci=="2"] <- "•"

pval2 <- p.adjust(rpval, "fdr")
pval2[rpval>0.1] <- "ns"
pval2[rpval<=0.1 & rpval>=0.05] <- "*"
pval2[rpval<=0.05 & rpval>=0.01] <- "**"
pval2[rpval<=0.01] <- "***"

# Create the barplot
g22 <- ggplot(data=rdfall2, aes(x=gnams, y=value, fill=variable)) +
  geom_bar(stat="identity",position = "stack", color = "#2A2A2A")+
  geom_text(aes(y=rep(rhfull,2), label=rep(paste0(round(rhval/rhfull,digits=2)*100,"%"),2)), hjust=0.3, vjust =1.3,
            color="black", size=3.5)+
  geom_text(aes(y=rep(rhfull,2), label=rep(pval2,2)), hjust=0.3,  vjust = -0.2,
            color="#2A2A2A", size=5)+
  geom_text(aes(y=rep(rep(0,length(rhfull)),2), label=rep(paste0("(",round(outd, digits = 2)*100, "%)"),2)), hjust=0.5,  vjust = -0.5,
            color="white", size=3.5)+
  scale_fill_manual(values=c("#998B84", "#D87141"), name="Targets shared\n with RAP-seq", aesthetics = "fill", labels=c('No', 'Yes'))+
  theme_classic2()+ rotate_x_text(45)+ylab("Number of targets in RAP-seq") + xlab("") +
  theme(text = element_text(size = 16),legend.position="top") 



#### PLOT 3
targ_gen <- "AQR"
meta_sub <- meta[gsub("-human","",meta$`Experiment target`) == targ_gen,]
file <- meta_sub$`File accession`[which(meta_sub$`File format`=="bigBed narrowPeak")]
countData_bb <- import.bb(paste0(file.path(paths$eclip_dir, ""),file,".bigBed"))
names(mcols(countData_bb)) <- c("name","score","signalValue","pValue","qValue","peak")

declip <- queryGff(queryRegions = countData_bb, gffData = annots) %>% as.data.table
#declip <- declip[which(declip$gene_type=="protein_coding"),]
#declip <- declip[which(declip$type=="gene"),]
#declip <- declip[which(as.numeric(p.adjust(10^(-as.numeric(declip$query_pValue)), method = "fdr")) <= 0.05),]


declip2 <- declip[which(declip$gene_name=="PES1"),]
coords <- lapply(lapply(strsplit(declip2$queryRange,":"), function(x) x[[2]]), function(x) strsplit(x, "-"))
coords <- do.call(rbind, lapply(coords, function(x) as.numeric(unlist(x))))
#  coords <- data.frame(start=declip2$start, end=declip2$end)

gr <- GRanges("chr22", IRanges(c(declip2$start[1], coords[,1]), c(declip2$start[1], coords[,2])), score=c(0, as.numeric(declip2$query_signalValue)))
tr <- new("track", dat=gr, type="data", format="BED")
setTrackStyleParam(tr, "color",  "#B63B37")
opts <- optimizeStyle(trackList(tr))

tr2 <- opts$tracks
#st2 <- opts$style
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                         org.Hs.eg.db,
                         gr=gr)


entrezID <- get("PES1", org.Hs.egSYMBOL2EG)
theTrack <- geneTrack(entrezID,TxDb.Hsapiens.UCSC.hg38.knownGene, asList = F)
setTrackStyleParam(theTrack, "color", "#3F5184")


trcks <- trackList(tr2, theTrack)


#setTrackStyleParam(trcks[[2]], name = "PES1")
#setTrackStyleParam(trcks[[1]], list(name = "eCLIP-seq"))
setTrackYaxisParam(trcks[[1]], "gp", list(cex=1))

setTrackStyleParam(trcks[[1]], "ylabgp", list(cex=1, col="#B63B37"))
setTrackStyleParam(trcks[[2]], "ylabgp", list(cex=1, col="#3F5184"))
setTrackXscaleParam(trcks[[1]], "draw", TRUE)

names(trcks) <- c("AQR\neCLIP score", "PES1")

vp <- viewTracks(trcks, chromosome="chr22", start=declip2$start[1], end=declip2$end[1]+1000, 
                 autoOptimizeStyle=TRUE) 
addGuideLine(c(30591851, 30591917), vp = vp, col = "#68769E", lwd = 2, lty = "dotted")
addGuideLine(c(30607047, 30607047), vp = vp, col = "#68769E", lwd = 2, lty = "dotted")
addArrowMark(list(x=30591917, 
                  y=2), # 2 means track 2 from the bottom.
             label="TSS (rfhg_229023.1,\nrfhg_229026.1)",
             col="#2A2A2A",
             vp=vp)

g3 <- recordPlot() 

##### plot 4

targ_gen <- "IGF2BP1"
file <- paste0(file.path(paths$rapseq_peaks_dir, ""),targ_gen, ".peaks.txt")
reclip <- read.table(file, header = T)
reclip <- reclip[which(reclip$gene_type=="protein_coding"),]
reclip2 <- reclip[which(reclip$gene_name=="PES1"),]
rcoords <- data.frame(reclip2$start, reclip2$end)

gr1 <- GRanges("chr22", IRanges(c(declip2$start[1], rcoords[,1]), c(declip2$start[1], rcoords[,2])), score=c(0, (as.numeric(reclip2$Rep1) + as.numeric(reclip2$Rep2))/2))
tr1 <- new("track", dat=gr1, type="data", format="BED")

setTrackStyleParam(tr1, "color",  "#3F844A")
opts1 <- optimizeStyle(trackList(tr1))
tr01 <- opts1$tracks
#trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene,
#                          org.Hs.eg.db,
#                         gr=gr1)


entrezID <- get("PES1", org.Hs.egSYMBOL2EG)
rtrack <- geneTrack(entrezID, TxDb.Hsapiens.UCSC.hg38.knownGene, asList = F)
setTrackStyleParam(rtrack, "color", "#3F5184")

rtrcks <- trackList(tr01, rtrack)


#setTrackStyleParam(trcks[[2]], name = "PES1")
#setTrackStyleParam(trcks[[1]], list(name = "eCLIP-seq"))
setTrackYaxisParam(rtrcks[[1]], "gp", list(cex=1))

setTrackStyleParam(rtrcks[[1]], "ylabgp", list(cex=1, col="#3F844A"))
setTrackStyleParam(rtrcks[[2]], "ylabgp", list(cex=1, col="#3F5184"))

setTrackXscaleParam(rtrcks[[1]], "draw", TRUE)

names(rtrcks) <- c("IGF2BP1\nRAPseq score", "PES1")

vp <- viewTracks(rtrcks, chromosome="chr22", start=declip2$start[1], end=declip2$end[1]+1000, 
                 autoOptimizeStyle=TRUE) 
addGuideLine(c(30591851,	30591917), vp = vp, col = "#68769E", lwd = 2, lty = "dotted")
addGuideLine(c(30607047, 30607047), vp = vp, col = "#68769E", lwd = 2, lty = "dotted")
addArrowMark(list(x=30591917, 
                  y=2), # 2 means track 2 from the bottom.
             label="TSS (rfhg_229023.1,\nrfhg_229026.1)",
             col="#2A2A2A",
             vp=vp)

g4 <- recordPlot() 

### together
alltrcks <- trackList(tr01, tr2,  rtrack)
setTrackYaxisParam(alltrcks[[1]], "gp", list(cex=1))
setTrackYaxisParam(alltrcks[[2]], "gp", list(cex=1))
setTrackStyleParam(alltrcks[[1]], "ylabgp", list(cex=1, col="#3F844A"))
setTrackStyleParam(alltrcks[[2]], "ylabgp", list(cex=1, col="#B63B37"))
setTrackStyleParam(alltrcks[[3]], "ylabgp", list(cex=1, col="#3F5184"))

setTrackXscaleParam(alltrcks[[1]], "draw", TRUE)
setTrackXscaleParam(alltrcks[[2]], "draw", TRUE)
names(alltrcks) <- c("IGF2BP1\nRAPseq score","AQR\neCLIP score",  "PES1")

vpall <- viewTracks(alltrcks, chromosome="chr22", start=declip2$start[1], end=declip2$end[1]+1000, 
                    autoOptimizeStyle=TRUE) 
addGuideLine(c(30591851,	30591917), vp = vp, col = "#68769E", lwd = 2, lty = "dotted")
addGuideLine(c(30607047, 30607047), vp = vp, col = "#68769E", lwd = 2, lty = "dotted")
addArrowMark(list(x=30591917, 
                  y=2), # 2 means track 2 from the bottom.
             label="TSS (rfhg_229023.1,\nrfhg_229026.1)",
             col="#2A2A2A",
             vp=vp)

g5 <- recordPlot() 

#################
ggsave(file=file.path(paths$output_dir, "figures/Fig1_net.svg"), plot=g1, width=12, height=9)
ggsave(file=file.path(paths$output_dir, "figures/Fig1_bar_eclip.svg"), plot=g2, width=5, height=10)
ggsave(file=file.path(paths$output_dir, "figures/Fig1_bar_rapseq.svg"), plot=g22, width=4, height=5)
svglite(file.path(paths$output_dir, "figures/Fig1_peak_aqr.svg"), width = 8, height = 4)
g3
dev.off()
svglite(file.path(paths$output_dir, "figures/Fig1_peak_igf2.svg"), width = 8, height = 4)
g4
dev.off()
svglite(file.path(paths$output_dir, "figures/Fig1_peak_eclip_rapseq.svg"), width = 8, height = 4)
g5
dev.off()
#################################################
#################################################
########   VALIDATION  ######## ######## ######## 
#################################################
#################################################
## libs #### libs #### libs ##
## libs #### libs #### libs ##


## reading files ##
## stage 1
gtex_pheno <- fread(paths$gtex_pheno_file)
lihc_surv <- fread(paths$lihc_survival_file)
genes_id <- fread(paths$gene_probe_map_file)
liver_ge2 <- readRDS(paths$liver_tcga_gtex_rds)
## stage 2
clsupp <- read_excel( paths$rbp_list_file)
cancerlitrbp <- read.table(paths$cancer_literature_rbp_file, sep = ',', header = TRUE)
cancerdegrbp <- read.table(paths$cancer_de_rbp_file, sep = ',', header = TRUE)
network_full <- read_tsv(paths$functional_network_file)
## stage 3
annots <- import.gff(paths$gencode_gtf_file)
meta <- read_tsv(paths$eclip_metadata_file)

## funs #### funs #### funs ##
discp <- function(p){
  pval <- p
  for(i in 1:length(pval)){
    if(is.nan(p[i])){
      pval[i] <- 0
    }else if(p[i]>0.1){
      pval[i] <- 0
    }else if(p[i]<=0.1 & p[i] > 0.05){
      pval[i] <- 1
    }else if(p[i]<=0.05 & p[i] > 0.01){
      pval[i] <- 2
    }else if(p[i]<=0.01 & p[i] > 0.001){
      pval[i] <- 3
    }else if(p[i]<=0.001){
      pval[i] <- 4
    }
  }
  return(pval)
}


discp2 <- function(p){
  pval <- p
  for(i in 1:length(pval)){
    if(is.na(p[i])){
      pval[i] <- "NA"
    }else if(p[i] == 0){
      pval[i] <- "ns"
    }else if(p[i] == 1){
      pval[i] <- "*"
    }else if(p[i] == 2){
      pval[i] <- "**"
    }else if(p[i] == 3){
      pval[i] <- "***"
    }else if(p[i] == 4){
      pval[i] <- "****"
    }
  }
  return(pval)
}

disccor <- function(c){
  c <- abs(round(c,digits=1))
  ccor <- c
  for(i in 1:length(ccor)){
    if(is.nan(c[i])){
      ccor[i] <- 0
    }else if(c[i] >= 0.7){
      ccor[i] <- 4
    }else if(c[i] < 0.7 &c[i] >= 0.5){
      ccor[i] <- 3
    }else if(c[i] < 0.5 & c[i] >= 0.3){
      ccor[i] <- 2
    }else if(c[i] < 0.3 & c[i] > 0){
      ccor[i] <- 1
    }
  }
  return(ccor)
}
## funs #### funs #### funs ##

################################################
############ READ YOUR DATA ####################
################################################
# read all edges and their weights
allEdgs <- readRDS(paths$all_interactions_rds)
colnames(allEdgs) <- c("fromi","toi","weight","from","to","color","arrows","arrows_type","k562")
allEdgs <- allEdgs[,c(4,5,3,9)]

allEdgs0 <- avgNets %>% as.matrix %>% t
allEdgs0[allEdgs0 < params$hepg2_network_threshold] <- 0 
rownames(allEdgs0) <- gsub(".1","",colnames(hepg2), fixed = T)[1:232]
colnames(allEdgs0) <- gsub(".1","",colnames(hepg2), fixed = T)[1:232]
allEdgs0 <- allEdgs0 %>% graph_from_adjacency_matrix(mode = "directed", weighted = TRUE)
Isolated <- which(igraph::degree(allEdgs0)==0)
allEdgs0 <- igraph::delete.vertices(allEdgs0, Isolated)
allEdgs <- data.frame(get.edgelist(allEdgs0), as.numeric(E(allEdgs0)$weight))
colnames(allEdgs) <- c("from", "to", "weight")

rownames(avgNetsk) <- gsub(".1","",colnames(k562), fixed = T)[1:232]
colnames(avgNetsk) <- gsub(".1","",colnames(k562), fixed = T)[1:232]

vk562 <- c()
for(i in 1:dim(allEdgs)[1]){
  nc <- which(geneNames_k562==allEdgs[i,]$from)
  nr <- which(geneNames_k562==allEdgs[i,]$to)
  if(length(nr)==0 | length(nc)==0){
    vk562[i] <- NA
  }else{
    vk562[i] <- avgNetsk[nr,nc]
  }
}

allEdgs$k562 <- vk562 

## survival analysis

g1tp <- c()
g2tp <- c()
psurv1 <- c()
psurv2 <- c()
ct1pval <- c()
ct2pval <- c()
cor1 <- c()
cor2 <- c()

allEdgs$from[allEdgs$from=="ATP5F1C"] <- "ATP5C1"
allEdgs$to[allEdgs$to=="ATP5F1C"] <- "ATP5C1"

allEdgs$from[allEdgs$from=="RACK1"] <- "GNB2L1"
allEdgs$to[allEdgs$to=="RACK1"] <- "GNB2L1"

# run survival analysis
for(i in 1:dim(allEdgs)[1]){
  genes <- c(allEdgs$from[i], allEdgs$to[i]) #c("RBM39","SRSF7")
  
  ################################################
  ################################################
  
  gene <- genes[1]
  gene_colid <- match(gene, colnames(liver_ge2))
  gene_ge <- liver_ge2[gene_colid]
  gene_ge <- as.data.frame(gene_ge)
  gene_ge$id <- substr(rownames(gene_ge),1,4)
  gene_ge$id[gene_ge$id=="TCGA"] <- "LIHC"
  gene_ge$id[gene_ge$id=="GTEX"] <- "control"
  if(length(which((table(gene_ge))>10))>=2){
    numToRem <- ((sort(table(gene_ge)))/sum((sort(table(gene_ge))))>=0.05)[((sort(table(gene_ge)))/sum((sort(table(gene_ge))))>=0.05)]
    gene_ge <- subset(gene_ge,!(as.numeric(unlist(gene_ge[1])) %in% as.numeric(names(numToRem))))
  }
  g1t <- t.test(gene_ge[gene_ge[,2]=="LIHC",1], gene_ge[gene_ge[,2]=="control",1])
  g1tp[i] <- g1t$p.value
  
  surv_gene <- lihc_surv[match(rownames(gene_ge), lihc_surv$sample) %>% na.omit] %>% as.data.frame()
  ge_d <- gene_ge[match(lihc_surv$sample, rownames(gene_ge)) %>% na.omit,1]
  ge_d2 <- rep("high",length(ge_d))
  ge_d2[ge_d <= median(ge_d)] <- "low"
  surv_gene$ge_status <- ge_d2
  
  #surv_gene <- surv_gene[which(surv_gene$OS.time<5*365),]
  fit <- survfit(Surv(OS.time, OS) ~ ge_status,
                 data = surv_gene)
  
  pval <- surv_pvalue(fit, surv_gene)
  if(pval$pval>0.1){
    pval_text <- "ns"
  }
  if(pval$pval<=0.1 & pval$pval > 0.05){
    pval_text <- "*"
  }
  if(pval$pval<=0.05 & pval$pval > 0.01){
    pval_text <- "**"
  }
  if(pval$pval<=0.01 & pval$pval > 0.001){
    pval_text <- "***"
  }
  if(pval$pval<=0.001){
    pval_text <- "****"
  }
  psurv1[i] <- pval$pval
  
  gene <- genes[2]
  gene_colid <- match(gene, colnames(liver_ge2))
  gene_ge <- liver_ge2[gene_colid]
  gene_ge <- as.data.frame(gene_ge)
  gene_ge$id <- substr(rownames(gene_ge),1,4)
  gene_ge$id[gene_ge$id=="TCGA"] <- "LIHC"
  gene_ge$id[gene_ge$id=="GTEX"] <- "control"
  if(length(which((table(gene_ge))>10))>=2){ # only when there are 2 sets of at least 10 the same values
    numToRem <- ((sort(table(gene_ge)))/sum((sort(table(gene_ge))))>=0.05)[((sort(table(gene_ge)))/sum((sort(table(gene_ge))))>=0.05)]
    gene_ge <- subset(gene_ge,!(as.numeric(unlist(gene_ge[1])) %in% as.numeric(names(numToRem))))
  }
  
  g2t <- t.test(gene_ge[gene_ge[,2]=="LIHC",1], gene_ge[gene_ge[,2]=="control",1])
  g2tp[i] <- g2t$p.value
  
  surv_gene <- lihc_surv[match(rownames(gene_ge), lihc_surv$sample) %>% na.omit] %>% as.data.frame()
  ge_d <- gene_ge[match(lihc_surv$sample, rownames(gene_ge)) %>% na.omit,1]
  ge_d2 <- rep("high",length(ge_d))
  ge_d2[ge_d <= median(ge_d)] <- "low"
  surv_gene$ge_status <- ge_d2
  
  #surv_gene <- surv_gene[which(surv_gene$OS.time<5*365),]
  fit <- survfit(Surv(OS.time, OS) ~ ge_status,
                 data = surv_gene)
  pval <- surv_pvalue(fit, surv_gene)
  if(pval$pval>0.1){
    pval_text <- "ns"
  }
  if(pval$pval<=0.1 & pval$pval > 0.05){
    pval_text <- "*"
  }
  if(pval$pval<=0.05 & pval$pval > 0.01){
    pval_text <- "**"
  }
  if(pval$pval<=0.01 & pval$pval > 0.001){
    pval_text <- "***"
  }
  if(pval$pval<=0.001){
    pval_text <- "****"
  }
  psurv2[i] <- pval$pval
  
  
  genes_colid <- match(genes, colnames(liver_ge2))
  genes_ge <- liver_ge2[genes_colid]
  
  colnames(genes_ge) <- c("gene1","gene2")
  if(length(which((table(genes_ge$gene1))>10))>=2 | length(which((table(genes_ge$gene2))>10))>=2){ # only when there are 2 sets of at least 10 the same values
    numToRem1 <- ((sort(table(genes_ge$gene1)))/sum((sort(table(genes_ge$gene1))))>=0.05)[((sort(table(genes_ge$gene1)))/sum((sort(table(genes_ge$gene1))))>=0.05)]
    numToRem2 <- ((sort(table(genes_ge$gene2)))/sum((sort(table(genes_ge$gene2))))>=0.05)[((sort(table(genes_ge$gene2)))/sum((sort(table(genes_ge$gene2))))>=0.05)]
    numToRem <- c(numToRem1, numToRem2)
    ntl1 <- !as.numeric(unlist(genes_ge$gene1)) %in% as.numeric(names(numToRem1))
    ntl2 <- !as.numeric(unlist(genes_ge$gene2)) %in% as.numeric(names(numToRem2))
    genes_ge <- subset(genes_ge, ntl1 & ntl2)
  }
  genes_ge$id <- substr(rownames(genes_ge),1,4)
  ct1 <- cor.test(genes_ge$gene1[genes_ge$id == "TCGA"], genes_ge$gene2[genes_ge$id == "TCGA"])
  ct1pval[i] <- ct1$p.value
  
  if(length(genes_ge$gene1[genes_ge$id == "GTEX"])==0){
    ct2pval[i] <- NaN
  }else{
    ct2 <- cor.test(genes_ge$gene1[genes_ge$id == "GTEX"], genes_ge$gene2[genes_ge$id == "GTEX"])
    ct2pval[i] <- ct2$p.value
  }
  
  cor1[i] <- ct1$estimate %>% as.numeric()
  cor2[i] <- ct2$estimate %>% as.numeric()
}
allEdgsSum <- allEdgs
allEdgsSum$deg1 <- g1tp
allEdgsSum$deg2 <- g2tp
allEdgsSum$survg1 <- psurv1
allEdgsSum$survg2 <- psurv2
allEdgsSum$corg1 <- cor1 # tcga
allEdgsSum$corg2 <- cor2 # gtex
allEdgsSum$corpg1 <- ct1pval # tcga
allEdgsSum$corpg2 <- ct2pval # gtex


### stage 2 ###
# adding 
edgs <- allEdgsSum

edgs$from[edgs$from=="ATP5C1"] <- "ATP5F1C"
edgs$to[edgs$to=="ATP5C1"] <- "ATP5F1C"

edgs$from[edgs$from=="GNB2L1"] <- "RACK1"
edgs$to[edgs$to=="GNB2L1"] <- "RACK1"

edgs$regulator_times_listed_as_rbp <- clsupp$Times_Listed_as_RBP[match(edgs$from, clsupp$gene_name)]
edgs$target_times_listed_as_rbp <- clsupp$Times_Listed_as_RBP[match(edgs$to, clsupp$gene_name)]
edgs$regulator_canonical <- clsupp$`canonical/non_canonical`[match(edgs$from, clsupp$gene_name)]
edgs$target_canonical <- clsupp$`canonical/non_canonical`[match(edgs$to, clsupp$gene_name)]
edgs$regulator_go_binding <- clsupp$Gene_Ontology_RNA_Binding[match(edgs$from, clsupp$gene_name)]
edgs$target_go_binding <- clsupp$Gene_Ontology_RNA_Binding[match(edgs$to, clsupp$gene_name)]
#literature rbps in cancer
edgs$regulator_cancer_lit <- cancerlitrbp$Number.of.Literatures[match(edgs$from, cancerlitrbp$Gene.Symbol)]
edgs$target_cancer_lit <- cancerlitrbp$Number.of.Literatures[match(edgs$to, cancerlitrbp$Gene.Symbol)]
#deg rbps in cancer
edgs$regulator_cancer_deg <- cancerdegrbp$Number.of.cancers[match(edgs$from, cancerdegrbp$Gene.Symbol)]
edgs$target_cancer_deg <- cancerdegrbp$Number.of.cancers[match(edgs$to, cancerdegrbp$Gene.Symbol)]

edgs$deg1 <- discp(edgs$deg1)
edgs$deg2 <- discp(edgs$deg2)
edgs$survg1 <- discp(edgs$survg1)
edgs$survg2 <- discp(edgs$survg2)
#edgs$corg1 <- disccor(edgs$corg1)
#edgs$corg2 <- disccor(edgs$corg2)
edgs$corpg1 <- discp(edgs$corpg1)
edgs$corpg2 <- discp(edgs$corpg2)

edgs$regulator_canonical <- as.numeric(as.factor(edgs$regulator_canonical)) # 1 is canonical
edgs$target_canonical <- as.numeric(as.factor(edgs$target_canonical)) # 1 is canonical

nams <- paste0(edgs$from,"-",edgs$to)
edgs_mat <- edgs[,-c(1:2)] 
rownames(edgs_mat) <- nams

edgs_mat2 <- edgs_mat %>% as.matrix
#heatmap(edgs_mat2 %>% scale, na_col = "grey90")

#edgs_mat2 <- edgs_mat[,1:6]
#edgs_mat2$cordiff <- abs(edgs_mat$corg1-edgs_mat$corg2)
#edgs_mat2$corpdiff <- abs(edgs_mat$corpg1-edgs_mat$corpg2)
#edgs_mat2$weight <- abs(edgs_mat2$weight*10)
#edgs_mat2$mns <- unname(rowMeans(edgs_mat2))
#pheatmap(edgs_mat2[edgs_mat2$mns>2,])

genes_ensm <- c()

for(i in 1:length(edgs$from)){
  genes <- c(edgs$from[i], edgs$to[i]) 
  symbols <- mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column="ENSEMBL")
  genes_ensm[i] <- paste(unname(symbols),sep="-",collapse = "-")
}


## FUNCOUP

genes_fc <- paste0(network_full$`2:Gene1`,"-",network_full$`3:Gene2`)
genes_fcrev <- paste0(network_full$`3:Gene2`,"-",network_full$`2:Gene1`)

m1 <- match(genes_fc, genes_ensm)
m1i <- which(!is.na(m1)==TRUE)

m2 <- match(genes_fcrev, genes_ensm)
m2i <- which(!is.na(m2)==TRUE)

m1m2i <- sort(c(m1i,m2i)) %>% unique

network_hepg2 <- network_full[m1m2i,]
# selecting subnetwork from FunCoup
symbols_1 <- mapIds(org.Hs.eg.db, keys = network_hepg2$`2:Gene1`, keytype = "ENSEMBL", column="SYMBOL")
symbols_2 <- mapIds(org.Hs.eg.db, keys = network_hepg2$`3:Gene2`, keytype = "ENSEMBL", column="SYMBOL")
network_hepg2$`2:Gene1` <- symbols_1 %>% unname
network_hepg2$`3:Gene2` <- symbols_2 %>% unname

edgs$from[edgs$from=="ATP5F1C"] <- "ATP5C1"
edgs$to[edgs$to=="ATP5F1C"] <- "ATP5C1"

edgs$from[edgs$from=="RACK1"] <- "GNB2L1"
edgs$to[edgs$to=="RACK1"] <- "GNB2L1"


fcogs <- paste0(network_hepg2$`2:Gene1`,"-",network_hepg2$`3:Gene2`)
fcogsrev <- paste0(network_hepg2$`3:Gene2`,"-",network_hepg2$`2:Gene1`)
ind <- list()
for(i in 1:dim(edgs)[1]){
  ogs <- paste0(edgs$from[i],"-",edgs$to[i])
  ind[[i]] <- c(which(fcogs == ogs), which(fcogsrev == ogs))
}
ind2 <- lapply(ind, function(x) length(x)) %>% unlist
ind2[which(ind2==1)] <- unlist(ind)
ind2[which(ind2==0)] <- NA

accs <- c()
#meths <- c("sd", "equal", "pretty", "quantile", "kmeans", "hclust", "bclust", "fisher", "jenks", "headtails")
#for(i in 1:length(meths)){
edgs$fc_score <- network_hepg2$`#0:PFC`[ind2]
edgs$weight <- abs(edgs$weight) # getting rid of the sign
#training
train_edgs <- edgs[which(!is.na(edgs$fc_score)),-c(1,2)]
test_edgs <- edgs[which(is.na(edgs$fc_score)),-c(1,2)]
rownames(train_edgs) <- paste0(edgs$from, '-', edgs$to)[which(!is.na(edgs$fc_score))]
rownames(test_edgs) <- paste0(edgs$from, '-', edgs$to)[which(is.na(edgs$fc_score))]


fc_score_binned <- classIntervals(train_edgs$fc_score, 2, style = "equal")
fc_score_binned$brks[1] <- 0
fc_score_binned2 <- cut(fc_score_binned$var, fc_score_binned$brks, labels=c("low","high"))


colnames(train_edgs) <- c("cfreq",
                          "k562i",
                          "regde",
                          "tarde",
                          "regsur",
                          "tarsurv",
                          "corLIHC",
                          "corGTEx",
                          "corpLIHC",
                          "corpGTEx",
                          "regtlrbp",
                          "tartlrbp",
                          "regcan",
                          "tarcan",
                          "reggob",
                          "targob",
                          "regclit",
                          "tarclit",
                          "regcdeg",
                          "tarcdeg",
                          "fc_score")

colnames(test_edgs) <- c("cfreq",
                         "k562i",
                         "regde",
                         "tarde",
                         "regsur",
                         "tarsurv",
                         "corLIHC",
                         "corGTEx",
                         "corpLIHC",
                         "corpGTEx",
                         "regtlrbp",
                         "tartlrbp",
                         "regcan",
                         "tarcan",
                         "reggob",
                         "targob",
                         "regclit",
                         "tarclit",
                         "regcdeg",
                         "tarcdeg",
                         "fc_score")
train_edgs$k562i <- as.matrix(train_edgs$k562i) %>% as.numeric
train_edgs$fc_score <- as.factor(fc_score_binned2)

## decision trees training
train_edgs$cfreq <- abs(train_edgs$cfreq)
# estimating FC score for unknown links with decision trees
model <- rpart(fc_score~., train_edgs, method = "class", minsplit = 5,
               minbucket = 2) # gini splitting, 
# display the results of modelling
printcp(model) 
plotcp(model)
summary(model)
rpart.plot(model, type=1, tweak = 1, Margin=0, fallen.leaves = F, extra=0) 

CV_pred <- predict(model, type="class")
true_pred <- train_edgs$fc_score %>% as.factor
names(true_pred) <- rownames(train_edgs)
cm <- confusionMatrix(CV_pred, true_pred) # CV accuracy 78%

## decision trees accuracy
unname(cm$overall["Accuracy"]) 

true_pred <- train_edgs$fc_score %>% as.factor
names(true_pred) <- rownames(train_edgs)

test_edgs$k562i <- as.matrix(test_edgs$k562i) %>% as.numeric
predict_model <- predict(model, test_edgs, type = "class") %>% as.data.frame
test_edgs$fc_score <- unname(as.matrix(predict_model))
#test_edgs_srt <- test_edgs[order(abs(test_edgs$cor_GTEx-test_edgs$cor_LIHC)*test_edgs$fc_score, decreasing = T),]

merge_traintest <- rbind(train_edgs, test_edgs)

train_edgs0 <- edgs[which(!is.na(edgs$fc_score)),-c(1,2)]
test_edgs0 <- edgs[which(is.na(edgs$fc_score)),-c(1,2)]
merge_traintest$fc_score_binned <- merge_traintest$fc_score
merge_traintest$fc_score <- c(train_edgs0$fc_score, test_edgs0$fc_score)

############################## stage 3 ############### eCLIP data

newrnms <- gsub("ATP5C1","ATP5F1C",rownames(merge_traintest))
newrnms <- gsub("GNB2L1","RACK1",newrnms)
rownames(merge_traintest) <- newrnms

# interactions for validation
inters <- merge_traintest
inters_names <- do.call(rbind, strsplit(rownames(inters),"-")) %>% as.data.frame
colnames(inters_names) <- c("regulator","target")

#allTargetGenes <- gsub("-human","",meta$`Experiment target`) %>% unique %>% sort

# reg to target
eclip_out <- c()

unr <- inters_names$regulator %>% unique
eclip_out <- rep(NA,length(inters_names$regulator))

for(i in 1:length(unr)){
  targ_gen <- unr[i]
  meta_sub <- meta[gsub("-human","",meta$`Experiment target`) == targ_gen,]
  
  if(isEmpty(meta_sub)){
    
  }else{
    
    file <- meta_sub$`File accession`[which(meta_sub$`File format`=="bigBed narrowPeak")]
    countData_bb <- import.bb(paste0(file.path(paths$eclip_dir, ""),file,".bigBed"))
    names(mcols(countData_bb)) <- c("name","score","signalValue","pValue","qValue","peak")
    
    overlaps <- as.data.table(queryGff(queryRegions = countData_bb, gffData = annots))
    overlaps <- overlaps[which(overlaps$gene_type=="protein_coding"),]
    overlaps <- overlaps[which(overlaps$type=="gene"),]
    overlaps <- overlaps[which(as.numeric(p.adjust(10^(-as.numeric(overlaps$query_pValue)), method = "fdr")) <= 0.05),]
    
    trgs <- inters_names$target[which(inters_names$regulator == targ_gen)]
    tmpecl <- c()
    for(j in 1:length(trgs)){
      tmpecl[j] <- length(which(overlaps$gene_name == trgs[j]))
    }
    eclip_out[which(inters_names$regulator == targ_gen)] <- tmpecl
  }
  print(i)
}

# target to reg

eclip_out2 <- c()

unt <- inters_names$target %>% unique
eclip_out2 <- rep(NA,length(inters_names$target))

for(i in 1:length(unt)){
  targ_gen <- unt[i]
  meta_sub <- meta[gsub("-human","",meta$`Experiment target`) == targ_gen,]
  
  if(isEmpty(meta_sub)){
    
  }else{
    
    file <- meta_sub$`File accession`[which(meta_sub$`File format`=="bigBed narrowPeak")]
    countData_bb <- import.bb(paste0(file.path(paths$eclip_dir, ""),file,".bigBed"))
    names(mcols(countData_bb)) <- c("name","score","signalValue","pValue","qValue","peak")
    
    overlaps <- as.data.table(queryGff(queryRegions = countData_bb, gffData = annots))
    overlaps <- overlaps[which(overlaps$gene_type=="protein_coding"),]
    overlaps <- overlaps[which(overlaps$type=="gene"),]
    overlaps <- overlaps[which(as.numeric(p.adjust(10^(-as.numeric(overlaps$query_pValue)), method = "fdr")) <= 0.05),]
    
    trgs <- inters_names$regulator[which(inters_names$target == targ_gen)]
    tmpecl <- c()
    for(j in 1:length(trgs)){
      tmpecl[j] <- length(which(overlaps$gene_name == trgs[j]))
    }
    eclip_out2[which(inters_names$target == targ_gen)] <- tmpecl
  }
  print(i)
}

inters$eclip_reg <- eclip_out
inters$eclip_targ <- eclip_out2

################################# RAPseq
# interactions for validation RAPseq

#allTargetGenes <- gsub("-human","",meta$`Experiment target`) %>% unique %>% sort

# reg to target
rapseq_out <- c()

unr <- inters_names$regulator %>% unique
rapseq_out <- rep(NA,length(inters_names$regulator))

raps <- list.files(paths$rapseq_peaks_dir)
allrs2 <- gsub("\\..*","",as.character(as.matrix(raps)))
gnsr <- c(allrs2, "IGF2BP2")

for(i in 1:length(unr)){
  targ_gen <- unr[i]
  meta_sub <- gnsr[!is.na(match(gnsr, targ_gen))]
  
  if(isEmpty(meta_sub)){
    
  }else{
    
    file <- paste0(file.path(paths$rapseq_peaks_dir, ""),targ_gen, ".peaks.txt")
    overlaps <- read.table(file, header = T)
    overlaps <- overlaps[which(overlaps$gene_type=="protein_coding"),]
    
    trgs <- inters_names$target[which(inters_names$regulator == targ_gen)]
    tmpecl <- c()
    for(j in 1:length(trgs)){
      tmpecl[j] <- length(which(overlaps$gene_name == trgs[j]))
    }
    rapseq_out[which(inters_names$regulator == targ_gen)] <- tmpecl
  }
  print(i)
}

# target to reg

rapseq_out2 <- c()
unt <- inters_names$target %>% unique
rapseq_out2 <- rep(NA,length(inters_names$target))

for(i in 1:length(unt)){
  targ_gen <- unt[i]
  meta_sub <- gnsr[!is.na(match(gnsr, targ_gen))]
  
  if(isEmpty(meta_sub)){
    
  }else{
    if(targ_gen == "IGF2BP2"){
      file1 <- paste0(file.path(paths$rapseq_peaks_dir, ""),targ_gen, "a.peaks.txt")
      overlaps1 <- read.table(file1, header = T)
      overlaps1 <- overlaps1[which(overlaps1$gene_type=="protein_coding"),]
      file2 <- paste0(file.path(paths$rapseq_peaks_dir, ""),targ_gen, "b.peaks.txt")
      overlaps2 <- read.table(file2, header = T)
      overlaps2 <- overlaps2[which(overlaps2$gene_type=="protein_coding"),]
      overlaps <- rbind(overlaps1, overlaps2)
    }else{
      file <- paste0(file.path(paths$rapseq_peaks_dir, ""),targ_gen, ".peaks.txt")
      overlaps <- read.table(file, header = T)
      overlaps <- overlaps[which(overlaps$gene_type=="protein_coding"),]
    }
    
    trgs <- inters_names$regulator[which(inters_names$target == targ_gen)]
    tmpecl <- c()
    for(j in 1:length(trgs)){
      tmpecl[j] <- length(which(overlaps$gene_name == trgs[j]))
    }
    rapseq_out2[which(inters_names$target == targ_gen)] <- tmpecl
  }
  print(i)
}

inters$rapseq_reg <- rapseq_out
inters$rapseq_targ <- rapseq_out2



#saveRDS(inters, file.path(paths$output_dir, "inters_hepg2_eclip_rapseq_new.rds"))
inters <- readRDS(file.path(paths$output_dir, "inters_hepg2_eclip_rapseq_new.rds"))
noc <- 3 # number of clusters

MclustNA <- function(x, n){
  x2 <- x
  wnix <- which(!is.na(x2))
  wzix <- which(x2 == 0)
  outvec <- rep(NA, length(x2))
  mc <- Mclust(as.numeric(na.omit(x2[-wzix])), n)$classification
  outvec[setdiff(wnix,wzix)] <- mc
  outvec[wzix] <- 0
  return(outvec)
}

edgs <- inters
edgs$k562i[edgs$k562i==0] <- NaN
edgsr <- data.frame('GRN_frequency' = abs(edgs$cfreq)/max(abs(edgs$cfreq)), # frequency from GRNs
                    'in_k562'= edgs$k562i/max(edgs$k562i, na.rm = T),
                    'DEG_LIHC' = rowMeans(data.frame(edgs$regde,edgs$tarde),na.rm = T)/max(rowMeans(data.frame(edgs$regde,edgs$tarde),na.rm = T), na.rm = T), # degs
                    'alter_survival' = rowMeans(data.frame(edgs$regsur,edgs$tarsurv),na.rm = T)/max(rowMeans(data.frame(edgs$regsur,edgs$tarsurv),na.rm = T), na.rm = T), # survival impact
                    'coexpression_change' = Mclust(abs(edgs$corLIHC-edgs$corGTEx), noc)$classification/noc, # difference between correlation
                    'times_listed_as_RBP' = Mclust(rowMeans(data.frame(edgs$regtlrbp,edgs$tartlrbp),na.rm = T),noc)$classification/noc, # times mentioned as RBP
                    'in_literature_as_cancer' = MclustNA(rowMeans(data.frame(edgs$regclit,edgs$tarclit), na.rm = T),2)/2, # literature
                    'in_cancers_as_DEG' = MclustNA(rowMeans(data.frame(edgs$regcdeg,edgs$tarcdeg), na.rm = T),2)/2, # cancer
                    'RBD'= rowMeans(data.frame(edgs$regcan,edgs$tarcan), na.rm = T)/max(rowMeans(data.frame(edgs$regcan,edgs$tarcan), na.rm = T)), # cancer
                    'regulator_eCLIP' = MclustNA(edgs$eclip_reg,2)/2,
                    'target_eCLIP' = MclustNA(edgs$eclip_targ,2)/2,
                    'regulator_RAPseq' = MclustNA(edgs$rapseq_reg,2)/2,
                    'target_RAPseq' = MclustNA(edgs$rapseq_targ,2)/2,
                    'fc_score'=edgs$fc_score_binned)

#rownames(edgsr) <- lapply(lapply(strsplit(rownames(edgs),"-"),function(x) c(x[2],x[1])), function(x) paste0(x,collapse="-")) %>% unlist
rownames(edgsr) <- gsub("-","→",rownames(edgs))

#res <- apply(edgsr[,-length(edgsr)], 1, function(x) sum(na.omit(x))) %>% as.data.frame # sum
res <- apply(edgsr[,-length(edgsr)], 1, function(x) mean(na.omit(x))) %>% as.data.frame # mean

res$score <- edgs$fc_score_binned %>% as.character
res$score[which(is.na(edgs$fc_score))] <- paste0(res$score[which(is.na(edgs$fc_score))],"*")
colnames(res) <- c("rank_score","fc_score")
#res$rank_score <- res$rank_score/max(res$rank_score)
edgsr$fc_score <- res$fc_score
scrs <- edgsr[order(res$rank_score, decreasing = T),]
scrs0 <- edgs[order(res$rank_score, decreasing = T),]
rsc <- res$rank_score[order(res$rank_score, decreasing = T)]

scrs00 <- scrs0[,c(1:20, 23:26 ,21:22)]
#write.csv2(scrs00, file = file.path(paths$output_dir, "supplementary_tables/Tab.S3_validationtable_3plus.csv"))
scrs <- read.csv2(file.path(paths$output_dir, "supplementary_tables/Tab.S3_validationtable_3plus.csv"))

#scrs2 <- scrs
#inters2 <- as.data.frame(do.call(rbind, strsplit(rownames(scrs2),"→")))
#colnames(inters2) <- c("regulator","target")
#inters2$grn_freq <- (abs(edgs$cfreq)*8)[order(res$rank_score, decreasing = T)]

#inters2$fcs <- scrs2$fc_score
####
bg <- clsupp$gene_name
bge <- mapIds(org.Hs.eg.db, bg, 'ENTREZID', 'SYMBOL')

thrs <- 1:dim(scrs)[1]
lstop <- list() 
lstopa <- list()

lstop2 <- list() 
lstopa2 <- list() 

### incremental threshold ###
for(i in 1:length(thrs)){
  genes <- unique(unlist(strsplit(scrs[1:thrs[i],1],"-"))) #unique(unlist(strsplit(rownames(scrs[1:thrs[i],]),"→")))
  genese <- mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL')
  edo <- enrichDGN(genese %>% unname, universe = bge %>% unname)
  eres <- edo@result
  #eres <- eres[p.adjust(eres$pvalue, method = "fdr") <= 0.1,]
  
  if(dim(eres)[1]==0){
    lstop[[i]] <- 0
    lstopa[[i]] <- 0
  }else{
    liv <- grep("liver",eres$Description)
    liv2 <- grep("Liver",eres$Description)
    hepa <- grep("hepato",eres$Description)
    hepa2 <- grep("Hepato",eres$Description)
    #canc <- grep("Cancer",eres$Description)
    #canc2 <- grep("cancer",eres$Description)
    #tumr <- grep("tumor",eres$Description)
    #tumr2 <- grep("Tumor",eres$Description)
    lall <- c(liv, hepa, liv2, hepa2) %>% unique #, canc, canc2, tumr, tumr2
    lstop[[i]] <- eres$pvalue[lall]
    names(lstop[[i]]) <- eres$Description[lall]
    lstopa[[i]] <- p.adjust(eres$pvalue[lall], "fdr")
    names(lstopa[[i]]) <- eres$Description[lall]
  }
  print(i)
}
pval_mean <- function(lstopa) {
  pl <- c()
  for(i in 1:length(lstopa)){
    if(lstopa[[i]][1]==0 | length(lstopa[[i]])<=1){
      pl[i] <- 1
    }else{
      pl[i] <- p.hmp(lstopa[[i]], L = length(lstopa[[i]])) %>% unname
    }
  }
  return(pl)
}

### sliding window ###
#win_siz <- 20
#for(i in 1:(length(thrs)-win_siz)){
#   genes <- unique(unlist(strsplit(rownames(scrs[i:(i+win_siz),]),"→")))
#   genese <- mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL')
#   edo <- enrichDGN(genese %>% unname, universe = bge %>% unname)
#   eres <- edo@result
#eres <- eres[p.adjust(eres$pvalue, method = "fdr") <= 0.1,]

#   if(dim(eres)[1]==0){
#     lstop2[[i]] <- 0
#     lstopa2[[i]] <- 0
#   }else{
#     liv <- grep("liver",eres$Description)
#     liv2 <- grep("Liver",eres$Description)
#     hepa <- grep("hepato",eres$Description)
#     hepa2 <- grep("Hepato",eres$Description)
#     lall <- c(liv, hepa, liv2, hepa2) %>% unique
#     lstop2[[i]] <- eres$pvalue[lall]
#     names(lstop2[[i]]) <- eres$Description[lall]
#     lstopa2[[i]] <- eres$p.adjust[lall]
#     names(lstopa2[[i]]) <- eres$Description[lall]
#   }
#   print(i)
# }

pval_top <- pval_mean(lstopa)
# pval_wind <- pval_mean(lstopa2)

### top interactions ###
#ptp <- which(pval_top <= 0.05)[length(which(pval_top <= 0.05))] #which(pval_top <= 0.05)[which(diff(which(pval_top <= 0.05))!=1)[1]]
ptp <- which.min(pval_top)
data1 <- data.frame(p=-log10(pval_top), thrs=thrs[1:length(pval_top)])
p1 <- ggplot(data=data1, aes(x=thrs,y=p)) +
  geom_bar(stat="identity", fill="#7A7A7A")+
  theme_minimal()+ 
  #geom_hline(yintercept=0.05, linetype="dashed", color = "#FE3200")+
  #geom_vline(xintercept=ptp, linetype="dashed", color = "#DAFF21", size=2)+
  #geom_text(aes(x=ptp, label=as.character(ptp), y=0.01), colour="#8F1003", vjust = 1.5, hjust=1.2) +
  #geom_text(aes(y=0.05, label=as.character(0.05), x=ptp), colour="#FE3200",vjust = -0.3, hjust = -0.3) +
  #geom_text(aes(y=0.5, label=formatC(min(pval_top), format = "e"), x=ptp), colour="#DAFF21",vjust = 0) +
  theme(legend.position = "none",text = element_text(size=15))+
  xlab("number of links")+
  ylab(bquote(-log[10]("P value"))) + 
  stat_peaks(col = "#0F61AF", span = 70, geom = "text_s", ignore_threshold = 0.05,
             arrow = arrow(length = grid::unit(1.5, "pt")), size=6, point.padding = 0.7)

nints <- 1:119
genes <- unique(unlist(strsplit(rownames(scrs[nints,]),"→"))) 
genese <- mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL')
edo <- enrichDGN(genese %>% unname, universe = bge %>% unname)

#write.table(sort(names(genese)),file.path(paths$output_dir, "top76_genes.txt"), quote = F, row.names = F, col.names = F)
#fcrange <- colorRampPalette(c("#f35c87","#5CB8F3"))
#p3 <- barplot(edo, showCategory=30)

data1 <- data.frame(pvalue = edo@result$p.adjust[which(edo@result$p.adjust<=0.05)], term = edo@result$Description[which(edo@result$p.adjust<=0.05)], count = as.numeric(unlist(lapply(strsplit(edo@result$GeneRatio,"/"), function(x) x[[1]]))[which(edo@result$p.adjust<=0.05)]))
#data1$genes <- factor(data1$genes, levels = data1$genes %>% unique)
data1$term <- factor(data1$term, levels = data1$term %>% unique)
#data1$pvalue <- -log10(data1$pvalue)

p3 <-   ggplot(data1, aes(x = term, y = count)) +
  geom_bar(aes(fill = pvalue), stat="identity")+ 
  scale_fill_gradient(low = "#f35c87", high = "#5CB8F3", na.value = NA, name = "P value")+
  coord_flip()+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50), limits = rev)+
  #ylab(bquote(-log[10](genes)))+
  xlab("") + theme_classic()+ theme(text = element_text(size = 16))+ ggtitle("DisGeNET")

ggarrange(p1, p3, nrow = 2, labels="AUTO")

### triangle plot ###
nn <- 119
scrstop <- scrs[1:nn,]
scrs130 <- scrs0[1:nn,]
rsc27 <- rsc[1:nn]
dfps <- data.frame(psurv1, psurv2)[order(res$rank_score, decreasing = T),]
ipvls <- apply(dfps, 1, p.hmp) 
ipvls2 <- ipvls[1:119] %>% unname

#nov2 <- (log10((scrs130$regclit+scrs130$tarclit)/2))*(log10((scrs130$regcdeg+scrs130$tarcdeg)/2))
nov2 <- (scrs130$regclit+scrs130$tarclit)/2

#nov2 <- log(nov2/min(nov2, na.rm = T)+2)+5

dfw <- data.frame(words=gsub("-","→",rownames(scrs130)),
                  freq=scrs130$cfreq,
                  coexp=log10(abs(scrs130$corLIHC/scrs130$corGTEx)),
                  surv= (log10((scrs130$regsur+scrs130$tarsurv)/2)),
                  nov=(log10((scrs130$regclit+scrs130$tarclit)/2))*(log10((scrs130$regcdeg+scrs130$tarcdeg)/2)))

dfw0 <- data.frame(words=gsub("-","→",rownames(scrs130)),
                   freq=scrs130$cfreq,
                   coexp=log2(abs(scrs130$corLIHC/scrs130$corGTEx)),
                   surv= ipvls2,
                   nov=nov2)
#dfw$freq <- dfw$freq/max(dfw$freq)
#dfw$coexp <- dfw$coexp/max(dfw$coexp)
#dfw$surv <- dfw$surv/max(dfw$surv)
#dfw$nov <- dfw$nov/max(dfw$nov,na.rm = T)


ng <- 3
labs <- dfw0$words
labs[(ng+1):length(labs)] <- "" # top 3

labs[which(dfw0$freq >=0.7)] <- dfw0$words[which(dfw0$freq >=0.7)]
tmpcoexp <- dfw0$coexp
tmpcoexp[dfw0$surv==0] <- 0
labs[c(which.max(tmpcoexp), which.min(tmpcoexp))] <- dfw$words[c(which.max(tmpcoexp), which.min(tmpcoexp))]
labs[which(dfw0$surv %in% dfw0$surv[order(dfw0$surv, decreasing = F)][1:3])] <- dfw0$words[which(dfw0$surv %in% (dfw0$surv[order(dfw0$surv, decreasing = F)][1:3]))]
labs
labs[c(which.min(dfw0$nov), which.max(dfw0$nov))] <- dfw0$words[c(which.min(dfw0$nov), which.max(dfw0$nov))]
tt <- list(
  size = 12,
  color = toRGB("black"),
  fill = "white")

axis <- function(title) {
  list(
    title = title,
    titlefont = list(
      size = 20
    ),
    tickfont = list(
      size = 15
    ),
    tickcolor = 'rgba(0,0,0,0)',
    ticklen = 5
  )
}

scale_values <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))}


fcs <- scrs130$fc_score
fcsb <- abs(as.numeric(as.factor(scrs130$fc_score_binned))-1)
sym <- rep(0, length(fcs))
sym[which(is.na(fcs))] <- 23
fcs[which(is.na(fcs))] <- fcsb[which(is.na(fcs))]

palette <- grDevices::colorRampPalette(c('darkcyan',"blueviolet"))(10)
clrs <- cut(fcs, 10, labels = palette)

tep <- rep("bottom center",nn)
tep[which(labs=="EWSR1→PES1")] <- "center right"
tep[which(labs=="HNRNPF→PES1")] <- "center right"
tep[which(labs=="HNRNPK→CKAP4")] <- "center right"
tep[which(labs=="NSUN2→PES1")] <- "center right"
tep[which(labs=="MSI2→PES1")] <- "center left"
tep[which(labs=="HSPD1→AKAP8L")] <- "center left"


lcol <- rep("darkgray",nn)
lcol[which(labs!="")] <- "orangered"

lwi <- rep(1,nn)
lwi[which(labs!="")] <- 2
bcolrs <- labs
bcolrs[bcolrs != ""] <- "#C33800"
bcolrs[bcolrs == ""] <- "#5F5F5F"

bstrok <- labs
bstrok[bstrok != ""] <- 2
bstrok[bstrok == ""] <- 1
dfw2 <- dfw0[which(dfw0$freq>=0.5),]
ggplot(dfw0, aes(x=coexp, y=surv, size=freq, fill=nov)) +
  geom_point(alpha=0.5, shape=21, color=bcolrs, stroke=bstrok %>% as.numeric) +
  scale_size(range = c(5, 20), name="Consensus frequency") +
  scale_fill_gradient(low = "#22FFCD", high = "#D1134A", na.value = NA, name = "Cancer literature\nmentions") +
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "#5F5F5F", linewidth=1) +
  theme_classic() +
  theme(text = element_text(size=16)) +
  ylab("Overall survival p value") + 
  xlab("Log2FC co-expression") + 
  geom_text_repel(aes(label=labs), size = 5,show_guide = T, max.overlaps = Inf)


#plot_ly() %>%
#  add_trace(
#   type='scatterternary',
#    a=scale_values(sqrt(dfw$surv+1)),
#    b=scale_values(sqrt(dfw$coexp+1)),
#    c=scale_values(sqrt(dfw$nov+1)),
#   mode = "markers+text",
#   text = gsub("→","↷",labs),
#   textposition = tep,
#   textfont = list(size=18),
#   marker = list( 
#     symbol = sym,
#     color = clrs,
#     opacity = 0.5,
#      size = dfw$freq*50,
#     line = list(color=lcol,'width' = lwi)
#   )) %>% layout(
#     ternary = list(
#       aaxis = axis('survival score'),
#       baxis = axis('co-expression FC'),
#       caxis = axis('popularity score')
#     )
#   )


#gt <- recordPlot() 
#svglite(file.path(paths$output_dir, "figures/Fig2_triangle.svg"), width = 8, height = 8)
#gt
#dev.off()

#################################################
scrstop <- scrs[1:nn,]

#scrstop <-  scrstop[which(scrstop$fc_score=="low*" | scrstop$fc_score=="low"),]

efrom <- lapply(strsplit(rownames(scrstop),"→"), function(x) x[1]) %>% unlist
eto <- lapply(strsplit(rownames(scrstop),"→"), function(x) x[2]) %>% unlist
el <- data.frame(efrom, eto)
escrs <- rsc27#[which(scrstop$fc_score=="low*" | scrstop$fc_score=="low")]

hig <- graph_from_data_frame(el, directed = TRUE)

hexprs <- read.table(paths$hepg2_expression_file)
hexprs1 <- (hexprs[,1:232] + hexprs[,233:464])/2

V(hig)$Node_degree = igraph::degree(hig)
mexps <- apply(hexprs, 1, median, na.rm=T)
mexps <- mexps[match(names(V(hig)), names(mexps))] %>% unname
V(hig)$mfc = mexps
#colfunc <- colorRampPalette(c("#C19DEE", "#FBAA11"))

vxnms2 <- names(V(hig))
tfrbp <- rep("RBP",length(vxnms2))
tfrbp[which(!is.na(match(vxnms2, TFs)))] <- "RBP and TF"
V(hig)$Node_shape = tfrbp

## sign of a link
alledgs <- get.edgelist(hig)
segds <- list()
tmpsign <- c()
pcr <- c()
scr <- c()
for(i in 1:dim(alledgs)[1]){
  for(j in 1:length(dfll10)){
    tmpnet <- dfll10[[j]] %>% t
    rownames(tmpnet) <- geneNames_hepg2
    colnames(tmpnet) <- geneNames_hepg2
    tmpsign[j] <- tmpnet[which(rownames(tmpnet) == alledgs[i,][1]),
                         which(colnames(tmpnet) == alledgs[i,][2])]
  }
  segds[[i]] <- tmpsign
  
  tmpgA <- hexprs1[which(rownames(hexprs1) == alledgs[i,][1]),] %>% as.matrix %>% as.numeric
  tmpgB <- hexprs1[which(rownames(hexprs1) == alledgs[i,][2]),] %>% as.matrix %>% as.numeric
  pcr[i] <- cor(tmpgA, tmpgB, method = 'pearson')
  scr[i] <- cor(tmpgA, tmpgB, method = 'spearman') ## 
}
msedgs <- do.call(rbind,segds)
msedgs2 <- msedgs[,-which(namns %in% c("CART","neunetreg"))]
msedgs2 <- cbind(msedgs2, scr %>% catnet %>% as.data.frame)

tmpscr <- scr %>% catnet
sgn <- c()
for(i in 1:dim(alledgs)[1]){
  tmpt <- table(msedgs2[i,] %>% as.matrix)
  tmpt <- tmpt[-which(names(tmpt) == 0)]
  
  if(length(tmpt)==1){
    sgn[i] <- names(tmpt) %>% as.numeric
  }else{
    if(tmpt[which(names(tmpt)=="-1")] == tmpt[which(names(tmpt)=="1")]){
      sgn[i] <- tmpscr[i]
    }else if(tmpt[which(names(tmpt)=="-1")] < tmpt[which(names(tmpt)=="1")]){
      sgn[i] <- 1
    }else{
      sgn[i] <- -1
    }
  } 
}

sgn <- sgn %>% as.character %>% as.factor
levels(sgn) <- c("activation","inhibition")
E(hig)$sign <- sgn


## eclip
ogns <- rep("", length(scrstop$regulator_eCLIP))
ogns[scrstop$regulator_eCLIP == 0] <- "•"
ogns[scrstop$regulator_eCLIP > 0] <- "★"

## rapseq
rogns <- rep("", length(scrstop$regulator_RAPseq))
rogns[scrstop$regulator_RAPseq == 0] <- "○"
rogns[scrstop$regulator_RAPseq > 0] <- "■"

E(hig)$weight <- escrs

E(hig)$label <- paste0(ogns, rogns)



gs1 <- ggraph(hig, layout = 'fr', circular=F) +
  geom_edge_arc(aes(width = weight, color = sign, label = label), label_colour = "#3F3F3F", label_size = 6, 
                alpha = 0.6, strength = 0.3, angle_calc = 'along', label_dodge = unit(2.5, 'mm'), 
                arrow = arrow(length = unit(6, 'pt')), end_cap = circle(6, 'pt')) + 
  geom_node_point(aes(size = Node_degree, color = mfc, shape = Node_shape)) + 
  geom_node_text(aes(label = name, size = 2), repel = TRUE, point.padding = unit(0.2, "lines")) +
  theme_graph(base_size = 14)+ 
  scale_size_continuous(range = c(3, 6)) +
  scale_colour_gradient2(name="Fold-change",
                         low = "#AA29EA",
                         mid = "#AAAAAA",
                         high = "#FF7E0C")+ 
  scale_edge_width(range = c(0.75, 2.25), name="Val. score") +
  scale_edge_color_manual(values=c("#FF320D","#0D91FF"), name="Interaction") +
  scale_label_size(range = c(1, 1)) +
  labs(size="Node degree", shape="Node shape")

ggsave(file=file.path(paths$output_dir, "figures/Fig2_net3_new.svg"), plot=gs1, width=12, height=9)


### nets high vs low
#################################################

scrstop <- scrs[1:nn,]
scrstop <-  scrstop[which(scrstop$fc_score=="high*" | scrstop$fc_score=="high"),]

outf <- strsplit(rownames(scrstop), "→") %>% unlist %>% unique
write.table(outf, file.path(paths$output_dir, "supplementary_tables/fchigh_rbps.txt"), quote = F, row.names = F, col.names = F)


efrom <- lapply(strsplit(rownames(scrstop),"→"), function(x) x[1]) %>% unlist
eto <- lapply(strsplit(rownames(scrstop),"→"), function(x) x[2]) %>% unlist
el <- data.frame(efrom, eto)
escrs <- rsc27[which(scrstop$fc_score=="high*" | scrstop$fc_score=="high")]

hig <- graph_from_data_frame(el, directed = TRUE)

hexprs <- read.table(paths$hepg2_expression_file)
hexprs1 <- (hexprs[,1:232] + hexprs[,233:464])/2

V(hig)$Node_degree = igraph::degree(hig)
mexps <- apply(hexprs, 1, median, na.rm=T)
mexps <- mexps[match(names(V(hig)), names(mexps))] %>% unname
V(hig)$mfc = mexps
#colfunc <- colorRampPalette(c("#C19DEE", "#FBAA11"))

vxnms2 <- names(V(hig))
tfrbp <- rep("RBP",length(vxnms2))
tfrbp[which(!is.na(match(vxnms2, TFs)))] <- "RBP and TF"
V(hig)$Node_shape = tfrbp


sgn <- scrstop$fc_score %>% as.character %>% as.factor
E(hig)$sign <- sgn

## eclip
ogns <- rep("", length(scrstop$regulator_eCLIP))
ogns[scrstop$regulator_eCLIP == 0] <- "•"
ogns[scrstop$regulator_eCLIP > 0] <- "★"

## rapseq
rogns <- rep("", length(scrstop$regulator_RAPseq))
rogns[scrstop$regulator_RAPseq == 0] <- "○"
rogns[scrstop$regulator_RAPseq > 0] <- "■"

E(hig)$weight <- escrs
E(hig)$label <- paste0(ogns, rogns)

gs1_high <- ggraph(hig, layout = 'fr', circular=F) +
  geom_edge_arc(aes(width = weight, color = sign, label = label), label_colour = "#3F3F3F", label_size = 6, 
                alpha = 0.6, strength = 0.3, angle_calc = 'along', label_dodge = unit(2.5, 'mm'), 
                arrow = arrow(length = unit(6, 'pt')), end_cap = circle(6, 'pt')) + 
  geom_node_point(aes(size = Node_degree, color = mfc, shape = Node_shape)) + 
  geom_node_text(aes(label = name, size = 2), repel = TRUE, point.padding = unit(0.2, "lines")) +
  theme_graph(base_size = 14)+ 
  scale_size_continuous(range = c(3, 6)) +
  scale_colour_gradient2(name="Fold-change",
                         low = "#AA29EA",
                         mid = "#AAAAAA",
                         high = "#FF7E0C")+ 
  scale_edge_width(range = c(0.75, 2.25), name="Val. score") +
  scale_edge_color_manual(values=c("#A2168B","#857282"), name="FunCoup5 score") +
  scale_label_size(range = c(1, 1)) +
  labs(size="Node degree", shape="Node shape")


## low

scrstop <- scrs[1:nn,]
scrstop <-  scrstop[which(scrstop$fc_score=="low*" | scrstop$fc_score=="low"),]

outf <- strsplit(rownames(scrstop), "→") %>% unlist %>% unique
write.table(outf, file.path(paths$output_dir, "supplementary_tables/fclow_rbps.txt"), quote = F, row.names = F, col.names = F)

efrom <- lapply(strsplit(rownames(scrstop),"→"), function(x) x[1]) %>% unlist
eto <- lapply(strsplit(rownames(scrstop),"→"), function(x) x[2]) %>% unlist
el <- data.frame(efrom, eto)
escrs <- rsc27[which(scrstop$fc_score=="low*" | scrstop$fc_score=="low")]

hig <- graph_from_data_frame(el, directed = TRUE)

hexprs <- read.table(paths$hepg2_expression_file)
hexprs1 <- (hexprs[,1:232] + hexprs[,233:464])/2

V(hig)$Node_degree = igraph::degree(hig)
mexps <- apply(hexprs, 1, median, na.rm=T)
mexps <- mexps[match(names(V(hig)), names(mexps))] %>% unname
V(hig)$mfc = mexps
#colfunc <- colorRampPalette(c("#C19DEE", "#FBAA11"))

vxnms2 <- names(V(hig))
tfrbp <- rep("RBP",length(vxnms2))
tfrbp[which(!is.na(match(vxnms2, TFs)))] <- "RBP and TF"
V(hig)$Node_shape = tfrbp


sgn <- scrstop$fc_score %>% as.character %>% as.factor
E(hig)$sign <- sgn

## eclip
ogns <- rep("", length(scrstop$regulator_eCLIP))
ogns[scrstop$regulator_eCLIP == 0] <- "•"
ogns[scrstop$regulator_eCLIP > 0] <- "★"

## rapseq
rogns <- rep("", length(scrstop$regulator_RAPseq))
rogns[scrstop$regulator_RAPseq == 0] <- "○"
rogns[scrstop$regulator_RAPseq > 0] <- "■"

E(hig)$weight <- escrs
E(hig)$label <- paste0(ogns, rogns)

gs1_low <- ggraph(hig, layout = 'fr', circular=F) +
  geom_edge_arc(aes(width = weight, color = sign, label = label), label_colour = "#3F3F3F", label_size = 6, 
                alpha = 0.6, strength = 0.3, angle_calc = 'along', label_dodge = unit(2.5, 'mm'), 
                arrow = arrow(length = unit(6, 'pt')), end_cap = circle(6, 'pt')) + 
  geom_node_point(aes(size = Node_degree, color = mfc, shape = Node_shape)) + 
  geom_node_text(aes(label = name, size = 2), repel = TRUE, point.padding = unit(0.2, "lines")) +
  theme_graph(base_size = 14)+ 
  scale_size_continuous(range = c(3, 6)) +
  scale_colour_gradient2(name="Fold-change",
                         low = "#AA29EA",
                         mid = "#AAAAAA",
                         high = "#FF7E0C")+ 
  scale_edge_width(range = c(0.75, 2.25), name="Val. score") +
  scale_edge_color_manual(values=c("#32EA89","#8ba698"), name="FunCoup5 score") +
  scale_label_size(range = c(1, 1)) +
  labs(size="Node degree", shape="Node shape")

ggsave(file=file.path(paths$output_dir, "figures/Fig2_net3_low.svg"), plot=gs1_low, width=10, height=8)
ggsave(file=file.path(paths$output_dir, "figures/Fig2_net3_high.svg"), plot=gs1_high, width=10, height=8)

### IGF2BP1

file1 <- paste0(file.path(paths$rapseq_peaks_dir, "IGF2BP1.peaks.txt"))
rigf1 <- read.table(file1, header = T)
rigf1 <- rigf1[which(rigf1$gene_type=="protein_coding"),]

meta_sub <- meta[gsub("-human","",meta$`Experiment target`) == "IGF2BP1",]
file2 <- meta_sub$`File accession`[which(meta_sub$`File format`=="bigBed narrowPeak")]
countData_bb <- import.bb(paste0(file.path(paths$eclip_dir, ""),file2,".bigBed"))
names(mcols(countData_bb)) <- c("name","score","signalValue","pValue","qValue","peak")
eigf1 <- as.data.table(queryGff(queryRegions = countData_bb, gffData = annots))
eigf1 <- eigf1[which(eigf1$gene_type=="protein_coding"),]
eigf1 <- eigf1[which(eigf1$type=="gene"),]
eigf1 <- eigf1[which(as.numeric(p.adjust(10^(-as.numeric(eigf1$query_pValue)), method = "fdr")) <= 0.05),]



allgnsi <- c(unique(eigf1$gene_name), unique(rigf1$gene_name)) %>% unique %>% sort
igns <- intersect(clsupp$gene_name, allgnsi)
igns2 <- mapIds(org.Hs.eg.db, igns, 'ENSEMBL', 'SYMBOL')
go_enrich <- enrichGO(gene = igns2 %>% unname,
                      universe =  clsupp$GRCh38.p7_ensembl_ID,
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      pAdjustMethod = "fdr",
                      ont = "BP",
                      pvalueCutoff = 0.01, 
                      qvalueCutoff = 0.2)

pthr <- 0.01
data1 <- data.frame(pvalue = go_enrich@result$p.adjust[which(go_enrich@result$p.adjust<=pthr)],
                    term = go_enrich@result$Description[which(go_enrich@result$p.adjust<=pthr)],
                    count = as.numeric(unlist(lapply(strsplit(go_enrich@result$GeneRatio,"/"), function(x) x[[1]]))[which(go_enrich@result$p.adjust<=pthr)]))
data1$term <- factor(data1$term, levels = data1$term %>% unique)


p11 <- ggplot(data1[1:20,], aes(x = term, y = count)) +
  geom_bar(aes(fill = pvalue), stat="identity")+ 
  scale_fill_gradient(low = "#f35c87", high = "#5CB8F3", na.value = NA, name = "P value")+
  coord_flip()+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50), limits = rev)+
  #ylab(bquote(-log[10](genes)))+
  xlab("") + theme_classic()+ theme(text = element_text(size = 16))+ ggtitle("GO Biological Process")


igns3 <- mapIds(org.Hs.eg.db, igns, 'ENTREZID', 'SYMBOL')
univ2 <- mapIds(org.Hs.eg.db, clsupp$gene_name, 'ENTREZID', 'SYMBOL')
igns3 <- igns3[-which(is.na(igns3))]
univ2 <- univ2[-which(is.na(univ2))]
go_enrich <- enrichDGN(gene = igns3 %>% unname,
                       universe =  univ2 %>% unname,
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.01, 
                       qvalueCutoff = 0.2)

pthr <- 0.01
data1 <- data.frame(pvalue = go_enrich@result$p.adjust[which(go_enrich@result$p.adjust<=pthr)],
                    term = go_enrich@result$Description[which(go_enrich@result$p.adjust<=pthr)],
                    count = as.numeric(unlist(lapply(strsplit(go_enrich@result$GeneRatio,"/"), function(x) x[[1]]))[which(go_enrich@result$p.adjust<=pthr)]))
data1$term <- factor(data1$term, levels = data1$term %>% unique)


p12 <- ggplot(data1[1:20,], aes(x = term, y = count)) +
  geom_bar(aes(fill = pvalue), stat="identity")+ 
  scale_fill_gradient(low = "#f35c87", high = "#5CB8F3", na.value = NA, name = "P value")+
  coord_flip()+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50), limits = rev)+
  #ylab(bquote(-log[10](genes)))+
  xlab("") + theme_classic()+ theme(text = element_text(size = 16))+ ggtitle("DisGeNET")

## KEGG
keggmed <- read.gmt(paths$msigdb_gmt_file)

keggen <- enricher(gene = igns,
                   universe =  clsupp$gene_name,
                   TERM2GENE=keggmed,
                   pAdjustMethod = "fdr")

pthr <- 0.01
data1 <- data.frame(pvalue =  keggen@result$p.adjust[which( keggen@result$p.adjust<=pthr)],
                    term =  keggen@result$Description[which( keggen@result$p.adjust<=pthr)],
                    count = as.numeric(unlist(lapply(strsplit( keggen@result$GeneRatio,"/"), function(x) x[[1]]))[which( keggen@result$p.adjust<=pthr)]))

trms <- gsub("_"," ",data1$term) %>% str_to_title
trms <- str_remove(trms, '(\\w+\\s+){1}')

data1$term <- factor(trms, levels = trms %>% unique)

p13 <- ggplot(data1[1:30,], aes(x = term, y = count)) +
  geom_bar(aes(fill = pvalue), stat="identity")+ 
  scale_fill_gradient(low = "#f35c87", high = "#5CB8F3", na.value = NA, name = "P value")+
  coord_flip()+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 60), limits = rev)+
  #ylab(bquote(-log[10](genes)))+
  xlab("") + theme_classic()+ theme(text = element_text(size = 16))+ ggtitle("Canonical Pathways")

p123 <- ggarrange(p13 + theme(text = element_text(size = 18)),
                  #p12 + theme(legend.position="bottom"),
                  p12 + theme(text = element_text(size = 18)),
                  nrow = 2, labels="AUTO")

ggsave(file=file.path(paths$output_dir, "figures/Fig2_enrich_supp.svg"), plot=p123, width=10, height=14)

ggsave(file=file.path(paths$output_dir, "figures/Fig2_gobp.svg"), plot=p11, width=8, height=5)
# arch plot

remotes::install_github("gastonstat/arcdiagram")
cols <- colorRampPalette(c("dimgray","seagreen1","darkmagenta"))(length(abs(log10(ecl_top130_sgn$pval))/max(abs(log10(ecl_top130_sgn$pval)))))
g <- graph_from_data_frame(ecl_top130_sgn , directed = TRUE, vertices = NULL)
edgelist <- get.edgelist(g)
ordn <- c(apply(edgelist, 1, as.matrix, collapse="")) %>% unique
degrees <- degree(g)
ordnn <- match(ordn,names(degrees))
degrees[ordnn]

arcplot(edgelist, sorted = F, horizontal = T, ordering =  V(g),
        lwd.arcs = ceiling(abs(log10(ecl_top130_sgn$pval))/10), 
        col.arcs = cols[rank(abs(log10(ecl_top130_sgn$pval))/10)], cex.nodes=log(degrees[ordnn])+0.5, col.nodes="dimgray")
a1 <- recordPlot()
plot.new() ## clean up device

############# peaks plot


#### PLOT 3
targ_gen <- "IGF2BP1"
meta_sub <- meta[gsub("-human","",meta$`Experiment target`) == targ_gen,]
file <- meta_sub$`File accession`[which(meta_sub$`File format`=="bigBed narrowPeak")]
countData_bb <- import.bb(paste0(file.path(paths$eclip_dir, ""),file,".bigBed"))
names(mcols(countData_bb)) <- c("name","score","signalValue","pValue","qValue","peak")

declip <- queryGff(queryRegions = countData_bb, gffData = annots) %>% as.data.table
declip <- declip[which(declip$gene_type=="protein_coding"),]
declip <- declip[which(declip$type=="gene"),]
declip <- declip[which(as.numeric(p.adjust(10^(-as.numeric(declip$query_pValue)), method = "fdr")) <= 0.05),]


# track with eCLIPseq I
declip2 <- declip[which(declip$gene_name=="PCBP2"),] #CCAR1
coords <- lapply(lapply(strsplit(declip2$queryRange,":"), function(x) x[[2]]), function(x) strsplit(x, "-"))
coords <- do.call(rbind, lapply(coords, function(x) as.numeric(unlist(x))))
gr <- GRanges(declip2$seqnames %>% unique %>% as.character, IRanges(c(declip2$start[1], coords[,1]), c(declip2$start[1], coords[,2])), score=c(0, as.numeric(declip2$query_signalValue)))
tr <- new("track", dat=gr, type="data", format="BED")
setTrackStyleParam(tr, "color",  "#B63B37")
opts <- optimizeStyle(trackList(tr))
tr2 <- opts$tracks

# track with eCLIPseq II
targ_gen <- "IGF2BP1"
file <- paste0(file.path(paths$rapseq_peaks_dir, ""),targ_gen, ".peaks.txt")
reclip <- read.table(file, header = T)
reclip <- reclip[which(reclip$gene_type=="protein_coding"),]
reclip2 <- reclip[which(reclip$gene_name=="PCBP2"),]
rcoords <- data.frame(reclip2$start, reclip2$end)
gr1 <- GRanges(declip2$seqnames %>% unique %>% as.character, IRanges(c(declip2$start[1], rcoords[,1]), c(declip2$start[1], rcoords[,2])), score=c(0, (as.numeric(reclip2$Rep1) + as.numeric(reclip2$Rep2))/2))
tr1 <- new("track", dat=gr1, type="data", format="BED")
setTrackStyleParam(tr1, "color",  "#3F844A")
opts1 <- optimizeStyle(trackList(tr1))
tr01 <- opts1$tracks


# track with gene
entrezID <- get("PCBP2", org.Hs.egSYMBOL2EG)
theTrack <- geneTrack(entrezID,TxDb.Hsapiens.UCSC.hg38.knownGene, asList = F)
setTrackStyleParam(theTrack, "color", "#3F5184")


trcks <- trackList(tr2, tr01, theTrack)


setTrackYaxisParam(trcks[[1]], "gp", list(cex=1))
setTrackYaxisParam(trcks[[2]], "gp", list(cex=1))

setTrackXscaleParam(trcks[[1]], "draw", TRUE)
setTrackXscaleParam(trcks[[2]], "draw", TRUE)

setTrackStyleParam(trcks[[1]], "ylabgp", list(cex=1, col="#B63B37"))
setTrackStyleParam(trcks[[2]], "ylabgp", list(cex=1, col="#3F844A"))
setTrackStyleParam(trcks[[3]], "ylabgp", list(cex=1, col="#3F5184"))



names(trcks) <- c("IGF2BP1\neCLIP score","IGF2BP1\nRAPseq score", "PCBP2")

vp <- viewTracks(trcks, chromosome=declip2$seqnames %>% unique %>% as.character, start=53441741, end=declip2$end[1]+1000, 
                 autoOptimizeStyle=TRUE) 
addGuideLine(c(53452096, 53452110), vp = vp, col = "#68769E", lwd = 2, lty = "dotted")
addGuideLine(c(53452303, 53452359), vp = vp, col = "#68769E", lwd = 2, lty = "dotted")
addGuideLine(c(53455886, 53456056), vp = vp, col = "#68769E", lwd = 2, lty = "dotted")
addGuideLine(c(53468773, 53468791), vp = vp, col = "#68769E", lwd = 2, lty = "dotted")
addArrowMark(list(x=53456056, 
                  y=2), # 2 means track 2 from the bottom.
             label="TSS (rfhg_151234.1,\nrfhg_151235.1,\nrfhg_151248.1)",
             col="#2A2A2A",
             vp=vp)

addArrowMark(list(x=53468791, 
                  y=2), # 2 means track 2 from the bottom.
             label="TSS (rfhg_151260.1)",
             col="#2A2A2A",
             vp=vp)

g6 <- recordPlot() 
svglite(file.path(paths$output_dir, "figures/Fig2_peak_eclip_rapseq.svg"), width = 8, height = 4)
g6
dev.off()


#################################################
#### stage 5 ####

### visualize results
topn <- nints %>% sort
#topn <- c(1:68, 134:155)# set 1 and 2
edgsr <- as.data.frame(edgsr)
ord_edgsr <- edgsr[as.numeric(order(res$rank_score, decreasing = T)),][topn,]
ord_edgsr[ord_edgsr==0] <- NaN
ord_edgsr <- ord_edgsr[,c(order(colMeans(ord_edgsr[,-length(ord_edgsr)]),decreasing = T), length(ord_edgsr))]
ord_edgsr <- ord_edgsr[,c(order(apply(ord_edgsr[,-length(ord_edgsr)], 2, function(x) length(which(is.na(x))))),length(ord_edgsr))]

## labels
ik562v <- edgs_full$k562i*10
ik562v[which(is.na(ik562v))] <- ""
edgs_full <- edgs[match(rownames(ord_edgsr), rownames(edgsr)),]
edgsr_labels <- data.frame('GRN_frequency' = format(abs(edgs_full$cfreq)*10), # frequency from GRNs
                           'in_k562' =  ik562v,
                           'DEG_LIHC' = paste0(discp2(edgs_full$regde),"→",discp2(edgs_full$tarde)), # degs
                           'alter_survival' = paste0(discp2(edgs_full$regsur),"→",discp2(edgs_full$tarsurv)), # survival impact
                           'coexpression_change' = paste0("C:",round(edgs_full$corLIHC, digits = 2),",H:",round(edgs_full$corGTEx, digits = 2)), # difference between correlation
                           'times_listed_as_RBP' = gsub("NA","Ø",paste0(round(edgs_full$regtlrbp, digits = 2),"→",round(edgs_full$tartlrbp, digits = 2))), # times mentioned as RBP
                           'in_literature_as_cancer' = gsub("NA","Ø",paste0(round(edgs_full$regclit, digits = 2),"→",round(edgs_full$tarclit, digits = 2))), # literature
                           'in_cancers_as_DEG' = gsub("NA","Ø",paste0(round(edgs_full$regcdeg, digits = 2),"→",round(edgs_full$tarcdeg, digits = 2))), # cancer
                           'RBD'= gsub("NA","Ø",paste0(gsub("2","nca",gsub("1","ca",edgs_full$regcan)),"→",(gsub("2","nca",gsub("1","ca",edgs_full$tarcan))))), # cancer
                           'regulator_eCLIP' = format(round(edgs_full$eclip_reg,digits = 1)),
                           'target_eCLIP' = format(round(edgs_full$eclip_targ,digits = 1)),
                           'regulator_RAPseq' = format(round(edgs_full$rapseq_reg,digits = 1)),
                           'target_RAPseq' = format(round(edgs_full$rapseq_targ,digits = 1)),
                           'fc_score' = edgs_full$fc_score_binned)
rownames(edgsr_labels) <- gsub("-","→",rownames(edgs_full))

##
edgsr2 <- ord_edgsr[,-length(ord_edgsr)] %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
edgsr2 <- edgsr2[order(edgsr2$colname),]

edgsr_labels2 <- edgsr_labels[,-length(ord_edgsr)] %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
edgsr_labels2 <- edgsr_labels2[order(edgsr_labels2$colname),]

edgsr2$label <- edgsr_labels2$value
edgsr2$colname2 <- factor(edgsr2$colname, levels = colnames(ord_edgsr)) 
edgsr2$rowname <- factor(edgsr2$rowname, levels = rev(unique(edgsr2$rowname)))
edgsr2 <- edgsr2[order(edgsr2$colname2),]
clrs <- rep("white",length(edgsr2$value))
clrs[is.na(edgsr2$value)] <- "grey30"
clrs[edgsr2$label=="NaN"] <- "white"

g1 <- ggplot(edgsr2, aes(x = colname2, y = rowname, fill = value), options(repr.plot.width = 10,repr.plot.height = 10)) +
  geom_tile(colour = "#636363", width=1) + theme_minimal()+ 
  scale_fill_gradient(high="#ff6200", low="#2275EC", na.value="white")+xlab("validation features")+
  annotate(geom = "text", x = 13.6,
           label = rev(ord_edgsr$fc_score),
           y = 1:length(topn), hjust = 0, col = "grey30")+
  theme(panel.grid.major = element_blank(),
        legend.position="left",
        legend.title=element_blank(),
        legend.margin = margin(0,-100,0,0),
        #plot.margin = unit(c(2, 2, 2, 2), "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(margin = unit(c(0,0,0,1),"npc")))+
  coord_cartesian(xlim = c(0, 15), clip = "off")+ylab("")+
  geom_text(aes(label = label), color = clrs, size = 3)+
  scale_x_discrete(labels=c("target_eCLIP" = "target eCLIP",
                            "regulator_eCLIP" = "regulator eCLIP",
                            "target_RAPseq" = "target RAPseq",
                            "regulator_RAPseq" = "regulator RAPseq",
                            "in_k562" = "link in k562",
                            "DEG_LIHC" = "DEG LIHC",
                            "times_listed_as_RBP" = "times listed as RBP",
                            "GRN_frequency" = "GRN frequency",
                            "coexpression_change" = "coexpression change",
                            "in_literature_as_cancer" = "in literature as cancer",
                            "in_cancers_as_DEG" = "in cancers as DEG",
                            "alter_survival" = "alter survival"))+
  geom_text(aes(x = 13.6,
                label = "FunCoup5 score",
                y = 0, angle = 45, fontface = c("plain")), hjust = 1, vjust = 1, color = "#636363", size = 3)


#ord_edgsr2 <- ord_edgsr[,-length(ord_edgsr)]
#datab <- apply(ord_edgsr2, 1, function(x) sum(na.omit(x))) %>% as.data.frame
#datab <- apply(ord_edgsr2, 1, function(x) mean(na.omit(x))) %>% as.data.frame

datab <- data.frame(val=res$rank_score[order(res$rank_score, decreasing = T)][topn])
rownames(datab) <- rownames(res)[order(res$rank_score, decreasing = T)][topn]
#datab$val <- datab$val#/max(res$rank_score)
datab$inters <- rownames(datab)
datab$inters <- factor(datab$inters, levels = unique(datab$inters))
g2 <- ggplot(data=datab, aes(x=rev(val), y=inters)) +
  geom_bar(stat="identity", colour="white", fill = "#636363")+ theme_minimal()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.005,0.6,0.07,-0.3), "npc"))+
  ylab("")+xlab("average score")

g12 <- ggarrange(g1, g2, ncol=2, labels=c("A",""), widths = c(1.2, 0.2))
#g123 <- ggarrange(g12, g3, ncol=1, nrow=2, labels=c("","B"))
#g123 + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
ggsave(g12, file = file.path(paths$output_dir, "figures/ggplot_val_table.svg"), limitsize = T, width = 25, height = 20, dpi = 700, bg = 'white')


#################################################
####### displaying a network in cytoscape #######
#################################################

# filter out some connections
displayNet2 <- displayNet
# remove detected with less than half methods
# displayNet2[round(abs(displayNet2)%%1,digits = 2) <= 0.5 & !abs(displayNet2)%in%c(1:20)] <- 0 
# only links that were validated with K562
# displayNet2[round(abs(displayNet2), digits = 2) <= 5] <- 0 

displayNet3 <- inter_hk[displayNet2!=0]

# plot how many links are common and uncommon between HepG2 and K562
val3 <- table(displayNet3, useNA = "always")
val4 <- melt(val3)
ggplot(val4[c(-which(is.na(val4$displayNet3))),], aes(x="", y=value, fill=displayNet3)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+ 
  theme_void() +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette="Set1")

# create directed network
net0 <- graph_from_adjacency_matrix(
  displayNet2,
  mode = c("directed"),
  weighted = T)

# construct edges
edges0 <- igraph::as_long_data_frame(net0)
cnms <- colnames(edges0)
cnms[4] <- "fromg"
cnms[5] <- "tog"
colnames(edges0) <- cnms
edges0$color <- rep("dimgray",length(edges0$from))
ewbin <- edges0$weight
ewbin[ewbin < 0] <- -1
ewbin[ewbin > 0] <- 1
aShape <- as.factor(ewbin)
levels(aShape) <- c("bar","arrow")
edges0$arrows.from.enabled <- TRUE
edges0$arrows.from.type <- as.character(aShape)
#edges0$weight <- igraph::as_data_frame(net0, what = c("edges"))$weight

# create directed network
net00 <- graph_from_adjacency_matrix(
  displayNet0,
  mode = c("directed"),
  weighted = T)

# construct edges
edges00 <- igraph::as_long_data_frame(net00)

we <- zeros(1,dim(edges0)[1])
vec0 <- paste0(edges00$from,edges00$to)
vec1 <- paste0(edges0$from,edges0$to)
we[match(vec0,vec1)] <- edges00$weight

edges0$k562 <- t(we) %>% as.data.frame()
saveRDS(edges0, file.path(paths$output_dir, "all_inters_hepg2.rds"))
#construct nodes
nodes <- data.frame(id=unique(c(edges0$fromg,edges0$tog)), shape='circle')
nodes$label <- unique(c(edges0$fromg,edges0$tog))

## type of interaction
inter <- rep(1,length(edges0$weight))
inter[edges0$weight<0] <- -1
inter <- as.factor(inter)
levels(inter) <- c(-1,1)
levels(inter) <- c("inhibits","activates")
inter2 <- as.character(inter)

#group <- rep("TG", length(nodes$id))
#group[which(nodes$id %in% TFs)] <- "TF"
hubs <- sort(table(c(edges0$fromg,
                     edges0$tog)),decreasing=T)

barplot(hubs[1:30], main="top 30 regulators", xlab="genes", ylab="degree",las=2)

set.seed(1)
hubsDegree <- discretize(unname(hubs), breaks = 3, labels=c("low","medium","high"), method = "cluster")
names(hubsDegree) <- names(hubs)
high_hubs <- names(hubsDegree[which(hubsDegree == "high")])
regu_hubs <- sort(table(edges$from[edges$from %in% high_hubs]), decreasing = T) 
regu_hubs_disc <- discretize(unname(regu_hubs), breaks = 2, labels=c("low","high"), method = "cluster")
high_regu_hubs <- regu_hubs[regu_hubs_disc=="high"]
barplot(high_regu_hubs, main="high regulators", xlab="", ylab="degree",las=2,col = "#2A6296")
hubsDegree2 <- rep("low",length(nodes$id))
hubsDegree2[match(names(high_regu_hubs),nodes$id)] <- "high"
#################################
## investigate high regulators ##
#################################

gostres <- gost(
  names(high_regu_hubs),
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("fdr"),
  domain_scope = c("annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
gostres$result$term_name[gostres$result$source=="KEGG"]
gostres$result$term_name[gostres$result$source=="WP"]
gostres$result$term_name[gostres$result$source=="REAC"]

# BP
simMatrix <- calculateSimMatrix(gostres$result$term_id[gostres$result$source=="GO:BP"],
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(gostres$result$p_value[gostres$result$source=="GO:BP"]), gostres$result$term_id[gostres$result$source=="GO:BP"])
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

treemapPlot(reducedTerms)
# MF
simMatrix <- calculateSimMatrix(gostres$result$term_id[gostres$result$source=="GO:MF"],
                                orgdb="org.Hs.eg.db",
                                ont="MF",
                                method="Rel")
scores <- setNames(-log10(gostres$result$p_value[gostres$result$source=="GO:MF"]), gostres$result$term_id[gostres$result$source=="GO:MF"])
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

treemapPlot(reducedTerms)

##############################################
############## create a network ##############
##############################################

group <- hubsDegree2 %>% unname %>% as.character
nodes2 <- data.frame(id=nodes$id,
                     group=group,
                     stringsAsFactors=FALSE)

edges2 <- data.frame(source=edges$from,
                     target=edges$to,
                     interaction=inter2,  # optional
                     weight=round(edges$weight, digits = 2), # numeric
                     stringsAsFactors=FALSE)

createNetworkFromDataFrames(nodes2, edges2, title="hepg2", collection="encode")

#########################################################
############## make a style to the network ##############
#########################################################

style.name = "myStyle02"
defaults <- list(NODE_SHAPE="ellipse",
                 NODE_SIZE=40,
                 EDGE_TRANSPARENCY=120)#,
#NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('Node Label','id','p')
#nodeLabelsColor <- mapVisualProperty('Node Label Color','group','d',c("low","medium","high"), c("black","black","white"))
arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
createVisualStyle(style.name, defaults, list(nodeLabels,arrowShapes))

allWei <- sort(as.numeric(unique(round(edges$weight,digits = 2))))
allWei2 <- allWei[which(allWei!=0)]
allWeiCols <- allWei2
allWeiCols[allWeiCols<0] <- "#FF320D"
allWeiCols[allWeiCols>0] <- "#0D91FF"

allWeiEdges <- rep("DASH_DOT", length(allWei2))
allWeiEdges[allWei2 >= 5 & allWei2 < 10] <- "SOLID"
allWeiEdges[allWei2 <= -5 & allWei2 > -10] <- "SOLID"
allWeiEdges[allWei2 >= 10] <- "ZIGZAG"
allWeiEdges[allWei2 <= -10] <- "ZIGZAG"

allWeiWi <- rep(1, length(allWei2))
allWeiWi[round(abs(allWei2)%%1,digits = 2) == 0 | round(abs(allWei2)%%1,digits = 2) == 0.1] <- 5

allWeiWi[round(abs(allWei2)%%1,digits = 2) >= 0.5 & round(abs(allWei2)%%1,digits = 2) <= 0.7] <- 4
allWeiWi[round(abs(allWei2)%%1,digits = 2) >= 0.3 & round(abs(allWei2)%%1,digits = 2) <= 0.5] <- 3
allWeiWi[round(abs(allWei2)%%1,digits = 2) >= 0.1 & round(abs(allWei2)%%1,digits = 2) <= 0.3] <- 2

setNodeColorMapping('group',c('low','medium','high'),colors=c('#dbdbdb','#808080','#212121'), style.name=style.name, mapping.type='d')

setEdgeLineWidthMapping('weight', allWei2, allWeiWi*2, style.name=style.name)
setEdgeColorMapping('weight',allWei2, allWeiCols, style.name=style.name)
setEdgeLineStyleMapping('weight',allWei2 %>% as.numeric, allWeiEdges, style.name=style.name)
setNodeLabelColorMapping('group', c("low","high"), c("#000000","#FFFFFF"), mapping.type = "d", style.name = style.name)
matchArrowColorToEdge(T, style.name = style.name)
setNodeFontSizeDefault(8, style.name = style.name)
setVisualStyle(style.name)

#####################################################
# check sparsity of a network
# all
mean(rowSums(apply(displayNet2,1,as.logical)))

# without zeros
csms <- rowSums(apply(displayNet2,1,as.logical))
mean(csms[csms!=0])

#### display a single hub with its connections ####
high_regu_hubs[1:3]
#hubName <- c("IGF2BP1","IGF2BP2","IGF2BP3") #ILF3, FUS
hubName <- c("AKAP1") #"NPM1","AKAP1","RBM39" 
nEdgesHub <- unique(c(which(!is.na(match(edges$from,hubName))),which(!is.na(match(edges$to,hubName)))))
edgesHub <- edges[nEdgesHub,]
nodes3 <- data.frame(id=nodes$id[match(unique(c(edgesHub$from,edgesHub$to)),nodes$id)],
                     #color=nodes$color, # categorical strings
                     group=group[match(unique(c(edgesHub$from,edgesHub$to)),nodes$id)],
                     stringsAsFactors=FALSE)
edges3 <- data.frame(source=edgesHub$from,
                     target=edgesHub$to,
                     interaction=inter2[nEdgesHub],  # optional
                     weight=edgesHub$weight, # numeric
                     stringsAsFactors=FALSE)
createNetworkFromDataFrames(nodes3, edges3, title="hepg2_hub", collection="encode")

write.table(nodes3$id,file.path(paths$output_dir, "geneNames_IGF2.txt"), quote = F, row.names = F, col.names = F)
setVisualStyle(style.name)
####################################################
gostres2 <- gost(
  nodes3$id,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.2,
  correction_method = c("fdr"),
  domain_scope = c("annotated"),
  custom_bg = NULL,
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
gostres2$result$term_name[gostres2$result$source=="WP"]
gostres2$result$term_name[gostres2$result$source=="KEGG"]
gostres2$result$term_name[gostres2$result$source=="HP"]
gostres2$result$term_name[gostres2$result$source=="HPA"]
gostres2$result$term_name[which(grepl("telom",gostres2$result$term_name))]
dftrms <- gostres2$result

simMatrix2 <- calculateSimMatrix(gostres2$result$term_id[gostres2$result$source=="GO:MF"],
                                 orgdb="org.Hs.eg.db",
                                 ont="MF",
                                 method="Rel")
scores2 <- setNames(-log10(gostres2$result$p_value[gostres2$result$source=="GO:MF"]), gostres2$result$term_id[gostres2$result$source=="GO:MF"])
reducedTerms2 <- reduceSimMatrix(simMatrix2,
                                 scores2,
                                 threshold=0.7,
                                 orgdb="org.Hs.eg.db")

treemapPlot(reducedTerms2)

#### SOME OTHER STUFF
# looking for subgraphs
net.bg <- graph_from_adjacency_matrix(
  displayNet2 %>% t,
  mode = c("undirected"))
net.bg <- simplify(net.bg)
Isolated = which(igraph::degree(net.bg)==0)
net.bg = delete.vertices(net.bg, Isolated)
V(net.bg)$size <- 6
V(net.bg)$frame.color <- "darkgray"
V(net.bg)$label <- "" 
E(net.bg)$arrow.mode <- 0
l <- layout_with_mds(net.bg)
fc <- cluster_fast_greedy(net.bg)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colors <- c("firebrick1","olivedrab4","lightblue2","navy","darkviolet","cyan1","thistle2","yellow1","lemonchiffon", "magenta1","darkorange","chartreuse")
V(net.bg)$color <- colors[fc$membership]
plot(net.bg, layout=l)

dfc <- data.frame(gene=fc$names,group=fc$membership)
lstgs <- dfc$gene[dfc$group==2]