#' @import AnnotationHub
#' @importFrom  babelgene orthologs
gmtPathways <- function(gmt.file,convertToEns,convertHu2Mm,verbose)
{
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  
  if (convertToEns & !convertHu2Mm)
  {
    tsmessage("... Retrieving gene annotation from AnnotationHub()",verbose = verbose)
    ah <- AnnotationHub::AnnotationHub()
    ahDb <- AnnotationHub::query(ah,pattern = c("Homo sapiens","EnsDb"), ignore.case = TRUE)
    id <- tail(rownames(mcols(ahDb)),n=1)
    edb <- ah[[id]]
    # Extract gene-level information from database
    ens.map <- subset(genes(edb,return.type = "data.frame"),seq_name%in%c(as.character(1:22),"X","Y","MT") & !gene_biotype%in%"LRG_gene")
    tsmessage(".. Start converting human symbols to human ensamble id",verbose = verbose)
    g = as.character(unique(unlist(pathways)))
    pathways = lapply(pathways, function(x,y=ens.map){r=y$gene_id[y$gene_name%in%x];r=r[!is.na(r)];return(unique(r))})
    tsmessage("Done!",verbose = verbose)
  }
  
  if (convertToEns & convertHu2Mm )
  {
    tsmessage(".. Start converting human symbols to mouse ensamble id",verbose = verbose)
    g = as.character(unique(unlist(pathways)))
    g = g[!startsWith(g,prefix = "ENSG")] 
    map.gene = babelgene::orthologs(genes = g,species = "mouse")
    pathways = lapply(pathways, function(x,y=map.gene) {r = unique(y$ensembl[y$human_symbol%in%x]);return(r[!is.na(r)])})
    tsmessage("Done!",verbose = verbose)
  }
  
  if (!convertToEns & convertHu2Mm )
  {
    tsmessage(".. Start converting human symbols to mouse symbols",verbose = verbose)
    g = as.character(unique(unlist(pathways)))
    g = g[!startsWith(g,prefix = "ENSG")]
    map.gene = babelgene::orthologs(genes = g,species = "mouse")
    pathways = lapply(pathways, function(x,y=map.gene) {r = unique(y$symbol[y$human_symbol%in%x]);return(r[!is.na(r)])})
    tsmessage("Done!",verbose = verbose)
  }
  
  pathways = pathways[sapply(pathways, length)>0]  
  pathways
}

#' Gene Set Enrichement Analysi on GF-ICF
#'
#' Compute GSEA for each cluster across a set of input pathways.
#' 
#' @param data list; GFICF object
#' @param gmt.file characters; Path to gmt file from MSigDB
#' @param nsim integer; number of simulation used to compute ES significance.
#' @param convertToEns boolean: Convert gene sets from gene symbols to Ensable id.
#' @param convertHu2Mm boolean: Convert gene sets from human symbols to Mouse Ensable id.
#' @param nt numeric; Number of cpu to use for the GSEA
#' @param minSize numeric; Minimal size of a gene set to test (default 15). All pathways below the threshold are excluded.
#' @param maxSize numeric; Maximal size of a gene set to test (default Inf). All pathways above the threshold are excluded.
#' @param verbose boolean; Show the progress bar.
#' @param seed integer; Seed to use for random number generation.
#' @param method string; Method to use GSEA or GSVA. Default is GSEA.
#' @return The updated gficf object.
#' @importFrom fgsea fgsea
#' @import fastmatch
#' @importFrom limma lmFit eBayes topTable
#' @import GSVA 
#' @export
runGSEA <- function(data,gmt.file,nsim=1000,convertToEns=T,convertHu2Mm=F,nt=2,minSize=15,maxSize=Inf,verbose=TRUE,seed=180582,method="GSEA")
{
  set.seed(seed)
  
  if (is.null(data$cluster.gene.rnk)) {stop("Please run clustcell function first")}
  mthod = base::match.arg(arg = method,choices = c("GSEA","GSVA"),several.ok = F)
  
  if (method == "GSEA")
  {
    tsmessage("Choosen method is GSEA...",verbose=verbose)
    data$gsea = list()
    data$gsea$pathways = gmtPathways(gmt.file,convertToEns,convertHu2Mm,verbose)
    data$gsea$es = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.gene.rnk))
    data$gsea$nes = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.gene.rnk))
    data$gsea$pval = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.gene.rnk))
    data$gsea$fdr = Matrix::Matrix(data = 0,nrow = length(data$gsea$pathways),ncol = ncol(data$cluster.gene.rnk))
    
    rownames(data$gsea$es) = rownames(data$gsea$nes) = rownames(data$gsea$pval) = rownames(data$gsea$fdr) = names(data$gsea$pathways)
    colnames(data$gsea$es) = colnames(data$gsea$nes) = colnames(data$gsea$pval) = colnames(data$gsea$fdr) = colnames(data$cluster.gene.rnk)
    
    pb = utils::txtProgressBar(min = 0, max = ncol(data$cluster.gene.rnk), initial = 0,style = 3)
    for (i in 1:ncol(data$cluster.gene.rnk))
    {
      df = as.data.frame(fgsea::fgseaMultilevel(pathways = data$gsea$pathways,stats = data$cluster.gene.rnk[,i],nPermSimple = nsim,gseaParam = 0,nproc = nt,minSize = minSize,maxSize = maxSize))[,1:7]
      data$gsea$es[df$pathway,i] = df$ES
      data$gsea$nes[df$pathway,i] = df$NES
      data$gsea$pval[df$pathway,i] = df$pval
      data$gsea$fdr[df$pathway,i] = df$padj
      utils::setTxtProgressBar(pb,i)
    }
    close(pb)
  
  data$gsea$stat = df[,c("pathway","size")]
  } else {
    tsmessage("Choosen method is GSVA...",verbose=verbose)
    data$gsva = list()
    data$gsva$pathways = gmtPathways(gmt.file,convertToEns,convertHu2Mm,verbose)
    data$gsva$DEpathways = NULL 
    data$gsva$res = Matrix::Matrix(data = 0,nrow = length(data$gsva$pathways),ncol = ncol(data$gficf))
    rownames(data$gsva$res) = names(data$gsva$pathways)
    colnames(data$gsva$res) = colnames(data$gficf)
    
    tsmessage("Start executiong GSVA cluster by cluster",verbose=verbose)
    options(warn=-1)
    u = unique(data$embedded$cluster)
    for (i in 1:length(u))
    {
      tsmessage(paste0("..Executing GSVA for cluster ",i," out of ",length(u)))
      cells = rownames(data$embedded)[data$embedded$cluster%in%u[i]]
      res = GSVA::gsva(expr = as.matrix(data$gficf[,cells]),gset.idx.list = data$gsva$pathways,kcdf="Gaussian",min.sz=minSize,max.sz=maxSize,parallel.sz=nt,method="gsva",verbose=F)
      data$gsva$res[rownames(res),cells] = res
      rm(res)
    }
    options(warn=0)
    data$gsva$res = data$gsva$res[armaRowSum(data$gsva$res!=0)>0,]
    
    tsmessage("Start executiong Limma cluster by cluster",verbose=verbose)
    for (i in 1:length(u))
    {
      tsmessage(paste0("..Calling DE pathways for cluster ",i," out of ",length(u)))
      clusters = data$embedded$cluster
      clusters[!clusters%in%u[i]] = "other"
      clusters[!clusters%in%"other"] = paste0("C",clusters[!clusters%in%"other"])
      design <- model.matrix(~ factor(clusters))
      colnames(design) <- c("ALL", paste0("C",u[i],"vsOTHER"))
      fit <- limma::lmFit(data$gsva$res, design)
      fit <- limma::eBayes(fit)
      df <- as.data.frame(limma::topTable(fit, coef=paste0("C",u[i],"vsOTHER"), number=Inf))
      df$pathway = rownames(df)
      df$cluster = u[i]
      data$gsva$DEpathways = rbind(data$gsva$DEpathways,df)
      rm(df)
    }
    rownames(data$gsva$DEpathways) = NULL
  }
  return(data)
}

#' Single cell Gene Set Enrichement Analysis on GF-ICF
#'
#' Compute GSEA for each cells across a set of input pathways by using NMF.
#' 
#' @param data list; GFICF object
#' @param geneID characters; The type of gene identifier to use, such as ensamble of symbol. 
#' @param species characters; Species name, such as human or mouse.
#' @param category characters; MSigDB collection abbreviation, such as H or C1.
#' @param subcategory characters; MSigDB sub-collection abbreviation, such as CGP or BP.
#' @param pathway.list list; Custom list of pathways. Each element correspond to a pathway a and contains a vector of genes.
#' @param nsim integer; number of simulation used to compute ES significance.
#' @param nt numeric; Number of cpu to use for the GSEA and NMF. Default is 0 (i.e., all available cores minus one)
#' @param minSize numeric; Minimal size of a gene set to test (default 15). All pathways below the threshold are excluded.
#' @param maxSize numeric; Maximal size of a gene set to test (default Inf). All pathways above the threshold are excluded.
#' @param verbose boolean; Show the progress bar.
#' @param seed integer; Seed to use for random number generation.
#' @param nmf.k numeric; Rank of NMF.
#' @param fdr.th numeric; FDR threshold for GSEA.
#' @param rescale string; If different by none, pathway's activity scores are resealed as Z-score. Possible values are none, byGS or byCell. Default is none.
#' @param normalization; normalization to use before to apply NMF. Possible values are gficf or cpm. Default and highly raccomanded is gficf.
#' @return The updated gficf object.
#' @importFrom fgsea fgsea
#' @import fastmatch
#' @importFrom RcppML nmf
#' @import utils
#' @import pointr
#' @import msigdbr
#' @importFrom BiocParallel SnowParam
#' @export
runScGSEA <- function(data,geneID,species,category,subcategory=NULL,pathway.list=NULL,nsim=10000,nt=0,minSize=15,maxSize=Inf,verbose=TRUE,seed=180582,nmf.k=100,fdr.th=0.05,gp=0,rescale="none",normalization="gficf")
{
  if(nt==0) {nt=detectCores()}
  #species = base::match.arg(arg = species,choices = c("human","mouse"),several.ok = F)
  geneID = base::match.arg(arg = geneID,choices = c("ensamble","symbol"),several.ok = F)
  rescale = base::match.arg(arg = rescale,choices = c("none","byGS","byCell"),several.ok = F)
  normalization = base::match.arg(arg = normalization,choices = c("gficf","cpm"),several.ok = F)
  options(RcppML.threads = nt)
  set.seed(seed)
  
  if (is.null(data$scgsea))
  {
    data$scgsea = list()
    if (normalization=="gficf")
    {
      if (!is.null(data$pca) && data$pca$type == "NMF"){
        if (data$dimPCA<nmf.k || data$pca$use.odgenes) {
          tsmessage("... Performing NMF",verbose=verbose)
          tmp = RcppML::nmf(data$gficf,k=nmf.k)
          data$scgsea$nmf.w <- Matrix::Matrix(data = tmp@w,sparse = T)
          data$scgsea$nmf.h <- t(Matrix::Matrix(data = tmp@h,sparse = T))
          rm(tmp);gc()
        } else {
          tsmessage(paste0("Found NMF reduction with k greaten or equal to ", nmf.k),verbose=T)
          pointr::ptr("tmp", "data$pca$genes")
          data$scgsea$nmf.w = tmp
          pointr::ptr("tmp2", "data$pca$cells")
          data$scgsea$nmf.h = tmp2
          rm(tmp,tmp2);gc()
        }
      } else {
        tsmessage("... Performing NMF",verbose=verbose)
        tmp = RcppML::nmf(data$gficf,k=nmf.k)
        data$scgsea$nmf.w <- Matrix::Matrix(data = tmp$w,sparse = T)
        data$scgsea$nmf.h <- t(Matrix::Matrix(data = tmp$h,sparse = T))
        rm(tmp);gc()
      }
    } else {
      tsmessage("... Performing NMF",verbose=verbose)
      tmp = RcppML::nmf(log1p(data$rawCounts),k=nmf.k)
      data$scgsea$nmf.w <- Matrix::Matrix(data = tmp$w,sparse = T)
      data$scgsea$nmf.h <- t(Matrix::Matrix(data = tmp$h,sparse = T))
      rm(tmp);gc()
    }
  } else {
    tsmessage("Found a previous scGSEA, thus the already computed NMF will be used",verbose=T)
    tsmessage("If you want to recompute NMF, please call resetScGSEA first",verbose=T)
    data$scgsea$es <- NULL
    data$scgsea$nes <- NULL
    data$scgsea$pval <- NULL
    data$scgsea$fdr <- NULL
    data$scgsea$pathways <- NULL
    data$scgsea$x <- NULL
    data$scgsea$stat <- NULL
  }
  
  tsmessage("Loading pathways...",verbose=verbose)
  if (!is.list(pathway.list)) {
    gs = msigdbr::msigdbr(species = species, category = category,subcategory = subcategory)
    if (geneID=="symbol") {
      data$scgsea$pathways = split(x = gs$gene_symbol, f = gs$gs_name)
    } else {
      data$scgsea$pathways = split(x = gs$ensembl_gene, f = gs$gs_name)
    }
  } else {
    data$scgsea$pathways = pathway.list; rm(pathway.list)
  }
  data$scgsea$es = Matrix::Matrix(data = 0,nrow = length(data$scgsea$pathways),ncol = ncol(data$scgsea$nmf.w))
  data$scgsea$nes = Matrix::Matrix(data = 0,nrow = length(data$scgsea$pathways),ncol = ncol(data$scgsea$nmf.w))
  data$scgsea$pval = Matrix::Matrix(data = 0,nrow = length(data$scgsea$pathways),ncol = ncol(data$scgsea$nmf.w))
  data$scgsea$fdr = Matrix::Matrix(data = 0,nrow = length(data$scgsea$pathways),ncol = ncol(data$scgsea$nmf.w))
    
  rownames(data$scgsea$es) = rownames(data$scgsea$nes) = rownames(data$scgsea$pval) = rownames(data$scgsea$fdr) = names(data$scgsea$pathways)
  
  tsmessage("Performing GSEA...",verbose=verbose)
  oldw <- getOption("warn")
  options(warn = -1)
  pb = utils::txtProgressBar(min = 0, max = ncol(data$scgsea$nmf.w), initial = 0,style = 3)
  
  # Optimize fgsea cpus assigning at least 100 pathways for cpus
  nt_fgsea <- ceiling(length(data$scgsea$pathways)/100)
  nt_fgsea <- ifelse(nt_fgsea>nt,nt,nt_fgsea)
  bpparameters <- BiocParallel::SnowParam(nt_fgsea)
  for (i in 1:ncol(data$scgsea$nmf.w))
  {
      df = as.data.frame(fgsea::fgseaMultilevel(pathways = data$scgsea$pathways,stats = data$scgsea$nmf.w[,i],nPermSimple = nsim,gseaParam = gp,BPPARAM = bpparameters,minSize = minSize,maxSize = maxSize))[,1:7]
      data$scgsea$es[df$pathway,i] = df$ES
      data$scgsea$nes[df$pathway,i] = df$NES
      data$scgsea$pval[df$pathway,i] = df$pval
      data$scgsea$fdr[df$pathway,i] = df$padj
      utils::setTxtProgressBar(pb,i)
  }
  base::close(pb)
  on.exit(options(warn = oldw))
  
  ix = is.na(data$scgsea$nes)
  if(sum(ix)>0) {
    data$scgsea$nes[ix] = 0
    data$scgsea$pval[ix] = 1
    data$scgsea$fdr[ix] = 1
  }
  
  data$scgsea$x = data$scgsea$nes
  data$scgsea$x[data$scgsea$x<0 | data$scgsea$fdr>=fdr.th] = 0
  data$scgsea$x = Matrix::Matrix(data = data$scgsea$nmf.h %*% t(data$scgsea$x),sparse = T)
  
  data$scgsea$stat = df[,c("pathway","size")]
  data$scgsea$x = data$scgsea$x[,armaColSum(data$scgsea$x)>0]
  
  if(rescale!="none"){
    if(rescale=="byGS") {
      data$scgsea$x = t(data$scgsea$x)
      data$scgsea$x = t( (data$scgsea$x - rowMeans(data$scgsea$x)) / apply(data$scgsea$x, 1, sd))
    }
    if(rescale=="byCell") {
      data$scgsea$x = (data$scgsea$x - rowMeans(data$scgsea$x)) / apply(data$scgsea$x, 1, sd)
    }
  }
  
  return(data)
}

#' Remove previous scGSEA analysis
#'
#' Remove previous scGSEA analysis
#' @param data list; GFICF object
#' @export
resetScGSEA <- function(data){
  data$scgsea <- NULL
  return(data)
}

