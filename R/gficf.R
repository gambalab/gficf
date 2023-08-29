#' Gene Frequency - Inverse Cell Frequency (GF-ICF)
#'
#' R implementation of the GF-ICF
#' Thanks to 3’-end scRNA-seq approaches, we can now have an accurate estimation of gene expression without having to account for gene length,
#' thus the number of transcripts (i.e. UMI) associated to each gene, strictly reflects the frequency of a gene in a cell, exactly like a word in a document.
#' GFICF (Gene Frequency - Inverse Cell Frequency) is analogous of TF-IDF scoring method as defined for tex mining With GFICF we consider a cell to be
#' analogous to a document, genes analogous to words and gene counts to be analogous of the word’s occurrence in a document.
#' 
#' @param M Matrix; UMI cell count matrix
#' @param QCdata list; QC cell object. 
#' @param cell_count_cutoff numeric; All genes detected in less than cell_count_cutoff cells will be excluded (default 5).
#' @param cell_percentage_cutoff2 numeric; All genes detected in at least this percentage of cells will be included (default 0.03, i.e. 3 percent of cells).
#' @param nonz_mean_cutoff numeric genes detected in the number of cells between the above mentioned cutoffs are selected only when their average expression in non-zero cells is above this cutoff (default 1.12).
#' @param  normalize logical; Rescale UMI counts before apply GFICF. Rescaling is done using EdgeR normalization.
#' @param storeRaw logical; Store UMI counts.
#' @param batches vector; Vector or factor for batch.
#' @param groups vector; Vector or factor for biological condition of interest.
#' @param filterGenes logical; Apply gene filter (default TRUE).
#' @param verbose boolean; Increase verbosity.
#' @param ... Additional arguments to pass to ComBat_seq call.
#' @return The updated gficf object.
#' 
#' @export
gficf = function(M=NULL,QCdata=NULL,cell_count_cutoff=5,cell_percentage_cutoff2=0.03,nonz_mean_cutoff=1.12,normalize=TRUE,storeRaw=TRUE,batches=NULL,groups=NULL,filterGenes=TRUE,verbose=TRUE,fastNomalization=F, ...)
{
  if(is.null(M) & is.null(QCdata)) {stop("Input data is missing!!")}
  
  data = list()
  if (!is.null(QCdata)) {
    data = QCdata
    rm(QCdata);gc(reset = T)
    if (!is.null(M)) {rm(M);gc()}
  } else {
    data$counts = M;rm(M);gc()
  }
  
  data = normCountsData(data,cell_count_cutoff,cell_percentage_cutoff2,nonz_mean_cutoff,normalize,batches,groups,verbose,filterGenes,fastNomalization, ...)
  data$gficf = tf(data$rawCounts,verbose = verbose)
  if (!storeRaw) {data$rawCounts=NULL;data$counts=NULL;gc()}
  data$w = getIdfW(data$gficf,verbose = verbose)
  data$gficf = idf(data$gficf,data$w,verbose = verbose)
  data$gficf = t(l.norm(t(data$gficf),norm = "l2",verbose = verbose))
  
  data$param <- list()
  data$param$cell_count_cutoff = cell_count_cutoff
  data$param$cell_percentage_cutoff2 = cell_percentage_cutoff2
  data$param$nonz_mean_cutoff = nonz_mean_cutoff
  data$param$normalized = normalize
  return(data)
}

#' @import Matrix
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom sva ComBat_seq
#' 
normCounts = function(M,cell_count_cutoff=5,cell_percentage_cutoff2=0.03,nonz_mean_cutoff=1.12,normalizeCounts=TRUE,batches=NULL,groups=NULL,verbose=TRUE,filterGene=FALSE,fastNomalization=FALSE, ...)
{
  if (filterGene) {
    tsmessage("Gene filtering..",verbose = verbose)
    M = filter_genes_cell2loc_style(data = M,cell_count_cutoff,cell_percentage_cutoff2,nonz_mean_cutoff)
  }
  
  if (normalizeCounts) 
  {
    if(!is.null(batches)){
      tsmessage("Correcting batches..",verbose = verbose)
      M = Matrix::Matrix(data = sva::ComBat_seq(counts = as.matrix(M),batch = batches,group = groups, ...),sparse = T)
    }
    tsmessage("Normalize counts..",verbose = verbose)
    if (!fastNomalization) {
      M <- Matrix::Matrix(edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts=M),normalized.lib.sizes = T)),sparse = T)
    } else {
      nt = ifelse(detectCores()>1,detectCores()-1,1)
      M <- scaleUMI(F,nt,F)
    }
  } else {
    data$rawCounts <- data$counts
    data$counts <- NULL; gc()
  }
  
  return(M)
}

#' @import Matrix
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom sva ComBat_seq
#' 
normCountsData = function(data,cell_count_cutoff=5,cell_percentage_cutoff2=0.03,nonz_mean_cutoff=1.12,normalizeCounts=TRUE,batches=NULL,groups=NULL,verbose=TRUE,filterGene=TRUE, fastNomalization=FALSE, ...)
{
  if (filterGene) {
    tsmessage("Gene filtering..",verbose = verbose)
    data$counts = filter_genes_cell2loc_style(data = data$counts,cell_count_cutoff,cell_percentage_cutoff2,nonz_mean_cutoff)
  }
  
  if (normalizeCounts) 
  {
    if(!is.null(batches)){
      tsmessage("Correcting batches..",verbose = verbose)
      data$counts = Matrix::Matrix(data = sva::ComBat_seq(counts = as.matrix(data$counts),batch = batches,group = groups, ...),sparse = T)
      gc()
    }
    tsmessage("Normalize counts..",verbose = verbose)
    if (!fastNomalization) {
      data$rawCounts <- Matrix::Matrix(edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts=data$counts),normalized.lib.sizes = T)),sparse = T)
    } else {
      nt = ifelse(detectCores()>1,detectCores()-1,1)
      data$rawCounts <- scaleUMI(data$counts,nt,F)
    }
  } else {
    data$rawCounts <- data$counts
    data$counts <- NULL; gc()
  }
  
  if(is.null(batches)){data$counts=NULL;gc()}
  
  return(data)
}

#' @import Matrix
#' 
tf = function(M,verbose)
{

  tsmessage("Apply GF transformation..",verbose = verbose)
  M =t(t(M) / armaColSum(M))
  
  return(M)
}

#' @import Matrix
#' 
idf = function(M,w,verbose)
{
  tsmessage("Applay ICF..",verbose = verbose)
  M = M[rownames(M) %in% names(w),]
  if(nrow(M)<length(w))
  {
    g = names(w)[!names(w)%in%rownames(M)]
    tmp = Matrix::Matrix(data = 0,nrow = length(g),ncol = ncol(M))
    rownames(tmp) = g
    colnames(tmp) = colnames(M)
    M = rbind(M,tmp)
  }
  M = M[names(w),]
  M = M * w
  return(M)
}

#' @import Matrix
#' 
getIdfW = function(M,type="classic",verbose)
{
  tsmessage("Compute ICF weigth..",verbose = verbose)
  nt = armaRowSum(M!=0)
  if (type == "classic") {w = log( (ncol(M)+1) / (nt+1) );rm(nt)}
  if (type == "prob") {w = log( (ncol(M) - nt) / nt );rm(nt)}
  if (type == "smooth") {w = log( 1 + ncol(M)/nt );rm(nt)}
  return(w)
}



l.norm = function (m, norm = c("l1", "l2"),verbose) 
{
  tsmessage(paste("Apply",norm),verbose = verbose)
  norm_vec = switch(norm, l1 = 1/armaRowSum(m), l2 = 1/sqrt(armaRowSum(m^2)))
  norm_vec[is.infinite(norm_vec)] = 0
  if (inherits(m, "sparseMatrix")) {
    norm_vec <- Diagonal(x = norm_vec)
    rownames(norm_vec) <- rownames(m)
    m <- norm_vec %*% m
  }
  else m * norm_vec
}

#' Save GFICF object.
#' 
#' Function to write a GFICF object to a file, preserving UMAP model.
#' 
#' @param data a GFICF object create by \code{\link{gficf}}.
#' @param file name of the file where the model is to be saved.
#' 
#' @examples 
#' # save
#' gficf_file <- tempfile("gficf_test")
#' saveGFICF(data, file = gficf_file)
#' 
#' # restore
#' data2 <- loadGFICF(file = gficf_file)
#' 
#' @export
saveGFICF <- function(data, file, verbose = TRUE) 
{
  wd <- getwd()
  tryCatch({
    # create directory to store files in
    mod_dir <- tempfile(pattern = "dir")
    #dir.create(mod_dir)
    gficf_dir <- file.path(mod_dir, "gficf")
    dir.create(gficf_dir,recursive = T)
    
    # save uwot object
    if (!is.null(data$reduction)) {
      if( data$reduction %in% c("umap","tumap"))
      {
        uwot_dir <- file.path(gficf_dir, "uwot")
        dir.create(uwot_dir,recursive = T)
        uwot_tmpfname <- file.path(uwot_dir,"uwot_obj")
        data$uwot <- uwot::save_uwot(model = data$uwot,file = uwot_tmpfname,verbose = verbose)
      }
    }
    
    # save gficf object
    gficfl_tmpfname <- file.path(gficf_dir, "data")
    saveRDS(data, file = gficfl_tmpfname)
    
    # archive the files under the temp dir into the single target file
    # change directory so the archive only contains one directory
    setwd(mod_dir)
    if(file.exists(file)) {file.rename(file,paste0(file,".bak"))}
    utils::tar(tarfile = file, files = "gficf/")
  },
  finally = {
    setwd(wd)
    if (file.exists(mod_dir)) {
      unlink(mod_dir, recursive = TRUE)
    }
    if(file.exists(paste0(file,".bak"))) {
      file.remove(paste0(file,".bak"))
    }
  })
}

#' Restore GFICF object.
#' 
#' Function to read a GFICF object from a file saved with \code{\link{saveGFICF}}.
#' 
#' @param file name of the file where the object is stored.
#' 
#' @examples 
#' gficf_file <- tempfile("gficf_data")
#' 
#' # restore
#' data2 <- loadGFICF(file = gficf_file)
#' 
#' @export
loadGFICF <- function(file,verbose = T)
{
  model <- NULL
  
  tryCatch({
    # create directory to store files in
    mod_dir <- tempfile(pattern = "dir")
    dir.create(mod_dir)
    utils::untar(file, exdir = mod_dir)
    
    # load gficf object
    gficf_fname <- file.path(mod_dir, "gficf/data")
    if (!file.exists(gficf_fname)) {
      stop("Can't find gficf data in ", file)
    }
    data <- readRDS(file = gficf_fname)
    
    # load umap obj if necessary
    if(data$reduction %in% c("umap","tumap"))
    {
      nn_fname <- file.path(mod_dir,"gficf/uwot/uwot_obj")
      if (!file.exists(nn_fname)) {
        stop("Can't find uwot object ", nn_fname, " in ", file)
      }
      data$uwot <- uwot::load_uwot(file = nn_fname, verbose = verbose)
    }
  },
  finally = {
    if (file.exists(mod_dir)) {
      unlink(mod_dir, recursive = TRUE)
    }
  })
  
  return(data)
}
