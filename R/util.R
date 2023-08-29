#' @import Matrix
#' 
scaleMatrix = function(x,rescale,centre)
{
  if (FALSE)
  {
    message("Rescaling..")
    bc_tot <- armaRowSum(x)
    median_tot <- stats::median(bc_tot)
    x <- base::sweep(x, 1, median_tot/bc_tot, '*')
    message("Rescaling Done!")
  }
  
  if (FALSE)
  {
    message("Centering data..")
    x <- base::sweep(x, 2, Matrix::colMeans(x), '-')
    x <- base::sweep(x, 2, base::apply(x, 2, sd), '/')
    message("Centering Done!")
  }
  return(x)
}


stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
tsmessage <- function(..., domain = NULL, appendLF = TRUE, verbose = TRUE,time_stamp = TRUE) {
  if (verbose) {
    msg <- ""
    if (time_stamp) {
      msg <- paste0(stime(), " ")
    }
    message(msg, ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}

#' Convert Ensamble IDs to Official Gene Symbols
#'
#' It uses biomart. If more the one gene is associated to the enamble, the first one retrived from
#' Biomart is used.
#' 
#' @param df data frame; Data frame containing the IDs to convert.
#' @param col characters; Name of column containing the ensamble ids.
#' @param organism characters; Organism of origin (i.e. human or mouse).
#' @param verbose boolean; Icrease verbosity.
#' @return The updated data frame with a new column called symb.
#'
#' @import AnnotationHub
#' @import fastmatch
#' 
#' @export
ensToSymbol = function(df,col,organism,verbose=T)
{
  organism = tolower(organism)
  organism = base::match.arg(arg = organism,choices = c("human","mouse"),several.ok = F)
  org.map = c("Homo Sapiens","Mus Musculus")
  names(org.map) = c("human","mouse")
  
  tsmessage("... Retrieving gene annotation from AnnotationHub()",verbose = verbose)
  ah <- AnnotationHub::AnnotationHub()
  ahDb <- AnnotationHub::query(ah,pattern = c(org.map[organism],"EnsDb"), ignore.case = TRUE)
  id <- tail(rownames(AnnotationHub::mcols(ahDb)),n=1)
  edb <- ah[[id]]
  ens.map <- subset(genes(edb,return.type = "data.frame"),seq_name%in%c(as.character(1:22),"X","Y","MT") & !gene_biotype%in%"LRG_gene")
  
  if(organism %in% "human")
  {
    tsmessage(".. Start converting human ensamble id to symbols",verbose = verbose)
    df$symb = NA
    df$symb = ens.map$gene_name[fastmatch::fmatch(df[,col],ens.map$gene_id)]
    tsmessage("Done!",verbose = verbose)
  }
  
  if (organism %in% "mouse")
  {
    tsmessage(".. Start converting mouse ensamble id to symbols",verbose = verbose)
    df$symb = NA
    df$symb = ens.map$gene_name[fastmatch::fmatch(df[,col],ens.map$gene_id)]
    tsmessage("Done!",verbose = verbose)
  }
  
  return(df)
}

#' Convert Official Gene Symbols to Ensamble IDs 
#'
#' It uses biomart. If more the one gene is associated to the enamble, the first one retrived from
#' Biomart is used.
#' 
#' @param df data frame; Data frame containing the IDs to convert.
#' @param col characters; Name of column containing the gene symbol.
#' @param organism characters; Organism of origin (i.e. human or mouse).
#' @param verbose boolean; Icrease verbosity.
#' @return The updated data frame with a new column called ens.
#'
#' @import AnnotationHub
#' @import fastmatch
#' 
#' @export
symbolToEns = function(df,col,organism,verbose=T)
{
  organism = tolower(organism)
  organism = base::match.arg(arg = organism,choices = c("human","mouse"),several.ok = F)
  org.map = c("Homo Sapiens","Mus Musculus")
  names(org.map) = c("human","mouse")
  
  tsmessage("... Retrieving gene annotation from AnnotationHub()",verbose = verbose)
  ah <- AnnotationHub::AnnotationHub()
  ahDb <- AnnotationHub::query(ah,pattern = c(org.map[organism],"EnsDb"), ignore.case = TRUE)
  id <- tail(rownames(AnnotationHub::mcols(ahDb)),n=1)
  edb <- ah[[id]]
  ens.map <- subset(genes(edb,return.type = "data.frame"),seq_name%in%c(as.character(1:22),"X","Y","MT") & !gene_biotype%in%"LRG_gene")
  
  if(organism %in% "human")
  {
    tsmessage(".. Start converting human symbols to human ensamble id",verbose = verbose)
    df$ens = NA
    df$ens = ens.map$gene_id[fastmatch::fmatch(df[,col],ens.map$symbol)]
    tsmessage("Done!",verbose = verbose)
  }
  
  if (organism %in% "mouse")
  {
    tsmessage(".. Start converting human symbols to mouse ensamble id",verbose = verbose)
    df$ens = NA
    df$ens = ens.map$gene_id[fastmatch::fmatch(df[,col],ens.map$symbol)]
    tsmessage("Done!",verbose = verbose)
  }
  
  return(df)
}


# Rcpp progress bar style
progress_for <- function(n, tot,display) {
  if (display) {
    message("0%   10   20   30   40   50   60   70   80   90   100%")
    message("[----|----|----|----|----|----|----|----|----|----|")
    # n:tot = nstars:50 -> nstars = (n*50)/tot
    nstars = floor((n*50)/tot)
    if(nstars>0)
      for (i in 1:nstars) {
        message("*", appendLF = FALSE)
        utils::flush.console()
      }
    message("|")
  }
}

# Filter cells with a cell2location style
# https://cell2location.readthedocs.io/en/latest/cell2location.utils.filtering.html
filter_genes_cell2loc_style = function(data,cell_count_cutoff=5,cell_percentage_cutoff2=0.03,nonz_mean_cutoff=1.12)
{
  nt = ifelse(detectCores()>1,detectCores()-1,1)
  data = data[armaRowSum(data)>0,]
  csums = armaRowSum(data!=0)
  gene_to_remove = csums <= cell_count_cutoff |  csums/ncol(data) <= cell_percentage_cutoff2
  gene_to_remove = gene_to_remove & armaColMeans(t(data),nt,FALSE)$mu1 <= nonz_mean_cutoff
  data = data[!gene_to_remove,]
  return(data)
}

#' Load in data from 10X
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#'
#' @param data.dir Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
#' files provided by 10X. A vector or named vector can be given in order to load
#' several data directories. If a named vector is given, the cell barcode names
#' will be prefixed with the name.
#' @param gene.column Specify which column of genes.tsv or features.tsv to use for gene names; default is 2
#' @param cell.column Specify which column of barcodes.tsv to use for cell names; default is 1
#' @param unique.features Make feature names unique (default TRUE)
#' @param strip.suffix Remove trailing "-1" if present in all cell barcodes.
#'
#' @return If features.csv indicates the data has multiple data types, a list
#'   containing a sparse matrix of the data from each type will be returned.
#'   Otherwise a sparse matrix containing the expression data will be returned.
#'
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # For output from CellRanger < 3.0
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' expression_matrix <- Read10X(data.dir = data_dir)
#' seurat_object = CreateSeuratObject(counts = expression_matrix)
#'
#' # For output from CellRanger >= 3.0 with multiple data types
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
#' data <- Read10X(data.dir = data_dir)
#' seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
#' seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)
#' }
#'
Read10X <- function(
    data.dir,
    gene.column = 1,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
) {
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, 'barcodes.tsv')
    gene.loc <- file.path(run, 'genes.tsv')
    features.loc <- file.path(run, 'features.tsv.gz')
    matrix.loc <- file.path(run, 'matrix.mtx')
    # Flag to indicate if this data is from CellRanger >= 3.0
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc) ) {
      stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    cell.barcodes <- read.table(file = barcode.loc, header = FALSE, sep = '\t', row.names = NULL)
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    } else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField,
        field = 1,
        delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    feature.names <- read.delim(
      file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning(
        'Some features names are NA. Replacing NA names with ID from the opposite column requested',
        call. = FALSE,
        immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[, gene.column]))
      replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column,
                    " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                    " Try setting the gene.column argument to a value <= to ", fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, gene.column])
    }
    # In cell ranger 3.0, a third column specifying the type of data was added
    # and we will return each type of data as a separate matrix
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) { # Return Gene Expression first
        lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
      }
      data <- lapply(
        X = lvls,
        FUN = function(l) {
          return(data[data_types == l, , drop = FALSE])
        }
      )
      names(x = data) <- lvls
    } else{
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    # Fix for Issue #913
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "CsparseMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}

detectCores <- function() {
  .Call("detectCoresCpp")
}

armaColSum <- function(M,nt=0,verbose=FALSE) {
  res = NULL
  c = class(M)
  if (nt==0) {nt = ifelse(detectCores()>1,detectCores()-1,1)}
  if(c[1]=="matrix") {
    res = armaColSumFull(M,nt,verbose)
  } else {
    if (c[1]=="lgCMatrix" || c[1]=="lgTMatrix") {
      res = Matrix::colSums(M)
    } else { 
      if(c[1]!="dgCMatrix") {M = as(M,"CsparseMatrix")}
      res = armaColSumSparse(M,nt,verbose)
    }
  }
  res = as.numeric(res)
  names(res) = colnames(M)
  return(res)
}

armaRowSum <- function(M,nt=0,verbose=FALSE) {
  return(armaColSum(t(M),nt,verbose))
}

#' @param M A sparse matrix from the Matrix package.
#' @param dir The directori in which to write the files. 
Write10X = function(M,dir) {
  if (!dir.exists(paths = dir)) {
    dir.create(path = dir,showWarnings = F,recursive = T)
  } 
  writeMMgz(x = M,file = paste0(dir,"/matrix.mtx.gz"))
  write.table(x = rownames(M),file = gzfile(paste0(dir,"/features.tsv.gz")),col.names = F,row.names = F)
  write.table(x = colnames(M),file = gzfile(paste0(dir,"/barcodes.tsv.gz")),col.names = F,row.names = F)
}

#' @param x A sparse matrix from the Matrix package.
#' @param file A filename that ends in ".gz".
writeMMgz <- function(x, file) {
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
      sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
    ),
    gzfile(file)
  )
  data.table::fwrite(
    x = summary(x),
    file = file,
    append = TRUE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
  )
}
