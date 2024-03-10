#' Create a Seurat Object
#'
#' Creating the Seurat object for bulk and single-cell ATAC-seq data.
#' Input takes a count matrix. The expected format of the input matrix is
#' features x cells. A set of genomic ranges (peaks or bins) must be supplied
#' along with the matrix, with the length of the ranges equal to the number of
#' rows in the matrix.
#'
#' @param matrix A count matrix.
#' @param ranges A set of GRanges corresponding to the rows of the input matrix.
#' @param meta.data A data.frame of meta data.
#' @param genome The name of a UCSC genome.
#' @param min.cells Include features detected in at least this many cells.
#'
#' @return A Seurat class object.
#'
#' @export

creatSeurat <- function(matrix, ranges, meta.data, genome, min.cells = 1){

  object <- CreateSeuratObject(
    counts = CreateChromatinAssay(counts = matrix,ranges = ranges,
                                  genome = genome,  min.cells = min.cells,
                                  min.features = 0),
    assay = "ATAC",min.cells = min.cells, min.features = 0,
    meta.data = meta.data
  )
  return(object)

}

#' Get the TF Motifs
#'
#' This function gets the TF Motifs from "JASPAR" or "cisBP" database.
#'
#' @param species The species source for the sequences, in Latin (Homo sapiens) or NCBI tax IDs (9606).
#' @param database A database of motif.
#'
#' @return A PFMatrixList object.
#'
#' @export

getMotifs <- function(species = "Homo sapiens",
                      database = c("JASPAR", "cisBP")[1], JASPAR = JASPAR2022){


  if(database == "JASPAR"){
    species <- "Homo sapiens"
    collection <- "CORE"
    opts <- list()
    opts["species"] <- species
    opts["collection"] <- collection
    JASPAR <- JASPAR
    out  <- TFBSTools::getMatrixSet(JASPAR, opts)

    if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
      names(out) <- paste(names(out), TFBSTools::name(out),
                          sep = "_")
    motifs <- out
  }
  if(database == "cisBP"){

    data("human_pwms_v1")
    sout <- sapply(strsplit(names(human_pwms_v1), split = "_"), function(s) c(s[3]))
    human_pwms_v2 <- human_pwms_v1[match(unique(sout), sout)]
    motifs <- human_pwms_v2
  }
  return(motifs)
}


#' PACells: identify clinical phenotype-associated cell states from Single-cell ATAC-seq Data
#'
#'  \code{PACells} is a novel approach to identify clinical phenotype-associated
#'  cell states from single-cell data using the phenotype, such as disease
#'  status, prognosis, and treatment response collected from bulk assays .
#'
#' @param sc_dataset A Seurat object of single-cell data.
#' @param bulk_dataset A Seurat object of bulk data.
#' @param phenotype A Phenotype annotation of bulk samples. It can be a binary
#' group indicator vector, continuous dependent variable, or clinical survival data.
#' @param motifs The TF Motifs from "JASPAR" or "cisBP" database.
#' @param cutoff Cutoff for the percentage of the PACells selected cells in total cells. This parameter is used to
#' restrict the number of the PACells selected cells. A smaller cutoff value (default \code{10\%}) is recommended
#' depending on the input data.
#' @param screenRatio Pre-screen candidate cells that are strongly related to the phenotype.
#' This parameter is used to pre-screen the percentage of candidate cells based Bcor measure.
#' An appropriate ratio (default \code{80\%}) is recommended depending on the input data.
#' @param family Response type for the regression model. It depends on the type of the given phenotype and
#' can be \code{family = gaussian} for linear regression, \code{family = binomial} for classification,
#' or \code{family = cox} for Cox regression.
#' @param group A vector of grouping information of single cells, which is calculated by the Signac package or provided by the user.
#' @param refgenome A reference genome (hg38 or hg19).
#'
#'
#' @return This function returns a Seurat-class object. It contains the results that cell state identified by PACells as strongly associated with the clinical phenotype.
#'
#' @export

PACells <- function(sc_dataset, bulk_dataset, phenotype, motifs,
                    cutoff = 0.1, screenRatio = 0.8,
                    family = c("binomial", "gaussian", "cox")[1],
                    group = NULL, sc_refgenome = c("hg38", "hg19")[1],
                    bulk_refgenome = c("hg38", "hg19")[1],seed = 123)
{


  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  #count matrix is the SummarizedExperiment class, motifs is the motif pwm matrix
  getTFscore <- function (object, motifs, refgenome = c("hg38", "hg19"))
  {
    if (refgenome == "hg19") {

      genom = BSgenome.Hsapiens.UCSC.hg19
    }
    if (refgenome == "hg38") {
      genom = BSgenome.Hsapiens.UCSC.hg38
    }
    count <- SummarizedExperiment(assays=SimpleList(counts = object@assays$ATAC@data),
                                  rowRanges = object@assays$ATAC@ranges, colData = object@meta.data)
    counts_GC <- addGCBias(count, genome = genom)
    counts_filtered <- filterPeaks(counts_GC)
    motif_ann <- matchMotifs(motifs, counts_filtered, genome = genom)
    dev <- computeDeviations(object = counts_filtered, annotations = motif_ann)
    z_score <- deviationScores(dev)
    return(z_score)
  }
  sc_TFscore <- getTFscore(sc_dataset, motifs, sc_refgenome)
  bulk_TFscore <- getTFscore(bulk_dataset, motifs, bulk_refgenome)
  cat("TF Activity Matrix done","\n")

  #get Similarty Matrix (KL)
  getSimilartyMatrix <- function(bulk_TFscore, sc_TFscore){

    dataset0 <- cbind(bulk_TFscore,sc_TFscore)
    dataset1 <- preprocessCore::normalize.quantiles(dataset0, keep.names = TRUE)

    bulk_TFscore1 <- dataset1[,1:dim(bulk_TFscore)[2]]
    sc_TFscore1 <- dataset1[,-c(1:dim(bulk_TFscore)[2])]

    nbins <- seq(min(dataset1), max(dataset1),length.out = dim(dataset1)[1]*2)###########
    sc_TFscore_count <- apply(sc_TFscore1, 2, function(x){hist(x, breaks = nbins, plot = FALSE)$counts})
    bulk_TFscore_count <- apply(bulk_TFscore1, 2, function(x){hist(x, breaks = nbins, plot = FALSE)$counts})

    bulk_TFscore_list <- lapply(seq_len(ncol(bulk_TFscore_count)), function(i) bulk_TFscore_count[,i])


    KlDist <- function(x, y, base=2){

      x[x == 0] <- 1e-15
      y[y == 0] <- 1e-15
      x <- x/sum(x)
      y <- y/sum(y)
      D1 <- sum(x * log(x/y, base = base))
      D2 <- sum(y * log(y/x, base = base))
      D <- (D1 + D2)/2
      return(list(D1 = D1, D2 = D2, D = D))

    }

    KLMatrix1 <- t(as.data.frame(lapply(bulk_TFscore_list,function(x){
      apply(sc_TFscore_count,2,function(y){KlDist(x,y)$D2})
    })))
    KL1Matrix <- max(KLMatrix1)-KLMatrix1
    rownames(KL1Matrix) <- colnames(bulk_TFscore)
    colnames(KL1Matrix) <- colnames(sc_TFscore)

    return(KL1Matrix)
  }
  simmtx <- getSimilartyMatrix(bulk_TFscore, sc_TFscore)
  cat("Similarty Matrix Done","\n")

  ####cell Pre-filtered
  ##X is all all single cell similarty matrix
  getCandidateCell <- function(simmtx, phenotype, screenRatio = 0.8,
                               family = "binomial"){
    X <- simmtx
    Y <- phenotype
    ratio <- screenRatio
    #    methods <- methods
    filter_num <-  dim(X)[2]*ratio
    if(family == "binomial"){
        fit <- sbisis(X = X, Y = Y,
                      candidate = filter_num,
                      method = "SBI-SIS-Pvalue",R=10)
    }
    if(family == "gaussian"){
        fit <- sbisis(X = X, Y = Y,
                      candidate = filter_num,
                      method = "SBI-SIS-Pvalue",R=10)
    }
    if(family == "cox"){
        fit <- sbisis.surv(X = X, Y = Y,
                           candidate = filter_num)
    }
    candidate_cell <- fit$candidate.set

    return(candidate_cell)
  }
  cand_cell <- getCandidateCell(simmtx, phenotype, screenRatio=screenRatio,
                                family = family)

  cat("Pre-screening Done","\n")

  ##X is all all single cell similarty matrix,
  ##Y: phenotype label,
  ##counts: the matrix of all cell
  ##sc_peak: the peaks of single cell data
  ##candidatecell: the candidate cells by runing the SBISIS
  ##ratio: The highest ratio of cells were identified
  ##resolution: the parameter findclusters()
  identifyCellSub <- function (X, Y, family = c("gaussian", "binomial", "cox")[2], object,
                               candidatecell, cutoff = 0.1, resolution = 1, group = NULL)
  {
    ratio <- cutoff
    seurat_object <- object
    seurat_object <- RunTFIDF(seurat_object, verbose = F)
    seurat_object <- FindTopFeatures(seurat_object, min.cutoff = "q10",
                                     verbose = F)
    seurat_object <- RunSVD(seurat_object, verbose = F)
    seurat_object <- RunUMAP(object = seurat_object, reduction = "lsi",
                             dims = 2:30, verbose = F)
    seurat_object <- FindNeighbors(object = seurat_object, reduction = "lsi",
                                   dims = 2:30, verbose = F)
    seurat_object <- FindClusters(object = seurat_object, verbose = FALSE,
                                  algorithm = 3, resolution = resolution)
    if(!is.null(group)){

      index <- group[candidatecell]###

    } else {
      index <- seurat_object$seurat_clusters[candidatecell]
    }

    inpute_x <- X[, candidatecell]

    if(family == "cox"){
      data = list(x = inpute_x, time = Y[,1], status = Y[,2])
      fit <- SGL::SGL(data, index, type = "cox", min.frac = 0.001,
                      nlam = 1000, standardize = F, alpha = 0, lambdas = seq(0.001,1,0.001))

    }
    if(family == "gaussian"){
      data = list(x= inpute_x, y = Y)
      fit <- SGL::SGL(data, index, type = "linear", min.frac = 0.001,
                      nlam = 1000, standardize = F, alpha = 0, lambdas = seq(0.001,1,0.001))

    }
    if(family == "binomial"){
      data = list(x= inpute_x, y = Y)
      fit <- SGL::SGL(data, index, type = "logit", min.frac = 0.001,
                      nlam = 1000, standardize = F, alpha = 0, lambdas = seq(0.001,1,0.001))

    }

    for (k in 1:dim(fit$beta)[2]) {
      nonzero_ratio <- (length(which((fit$beta[, k] > 0))))/(dim(X)[2])
      num <- k
      if (nonzero_ratio < ratio)
        break
    }

    sgl_results <- fit$beta[, num]
    label_sgl <- rep("Background", dim(X)[2])
    names(label_sgl) <- colnames(X)
    label_sgl[colnames(inpute_x)[which(sgl_results > 0)]] = "PACells"
    cat("Label Summary (B/P): ", table(label_sgl))
    sc_data <- seurat_object
    sc_data@meta.data$PACells_label <- factor(label_sgl, levels = c("PACells",
                                                                    "Background"))
    return(sc_data)
  }

  sc_data <-  identifyCellSub(X=simmtx, Y=phenotype, family = c("gaussian", "binomial", "cox")[2],
                              object = sc_dataset, candidatecell=cand_cell,
                              cutoff = cutoff, resolution = 1)

  return(sc_data)
}


#' PACells.RNA: identify clinical phenotype-associated cell states from Single-cell RNA-seq Data
#'
#'  \code{PACells} is a novel approach to identify clinical phenotype-associated
#'  cell states from single-cell data using the phenotype, such as disease
#'  status, prognosis, and treatment response collected from bulk assays .
#'
#' @param sc_dataset A expression matrix of single-cell data.
#' @param bulk_dataset A expression matrix of bulk data.
#' @param phenotype A Phenotype annotation of bulk samples. It can be a binary
#' group indicator vector, continuous dependent variable, or clinical survival data.
#' @param cutoff Cutoff for the percentage of the PACells selected cells in total cells. This parameter is used to
#' restrict the number of the PACells selected cells. A smaller cutoff value (default \code{10\%}) is recommended
#' depending on the input data.
#' @param screenRatio Pre-screen candidate cells that are strongly related to the phenotype.
#' This parameter is used to pre-screen the percentage of candidate cells based Bcor measure.
#' An appropriate ratio (default \code{80\%}) is recommended depending on the input data.
#' @param family Response type for the regression model. It depends on the type of the given phenotype and
#' can be \code{family = gaussian} for linear regression, \code{family = binomial} for classification,
#' or \code{family = cox} for Cox regression.
#' @param group A vector of grouping information of single cells, which is calculated by the Seurat package or provided by the user.
#'
#'
#' @return This function returns a Seurat-class object. It contains the results that cell state identified by PACells as strongly associated with the clinical phenotype.
#'
#' @export

PACells.RNA <- function(sc_dataset, bulk_dataset, phenotype,
                    cutoff = 0.1, screenRatio = 0.8,
                    family = c("binomial", "gaussian", "cox")[1],
                    group = NULL, seed = 123)
{

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  #get Similarty Matrix (KL)
  getSimilartyMatrix <- function(bulk_dataset, sc_dataset){
    common <- intersect(rownames(bulk_count), rownames(sc_count))


    dataset0 <- cbind(bulk_count[common,],sc_count[common,])
    dataset1 <- preprocessCore::normalize.quantiles(dataset0, keep.names = TRUE)

    bulk_norm <- dataset1[,1:dim(bulk_count)[2]]
    sc_norm <- dataset1[,-c(1:dim(bulk_count)[2])]

    nbins <- seq(min(dataset1), max(dataset1),length.out = dim(dataset1)[1]*2)###########
    sc_1 <- apply(bulk_norm, 2, function(x){hist(x, breaks = nbins, plot = FALSE)$counts})
    bulk_1 <- apply(bulk_norm, 2, function(x){hist(x, breaks = nbins, plot = FALSE)$counts})

    bulk_1_list <- lapply(seq_len(ncol(bulk_1)), function(i) bulk_1[,i])


    KlDist <- function(x, y, base=2){

      x[x == 0] <- 1e-15
      y[y == 0] <- 1e-15
      x <- x/sum(x)
      y <- y/sum(y)
      D1 <- sum(x * log(x/y, base = base))
      D2 <- sum(y * log(y/x, base = base))
      D <- (D1 + D2)/2
      return(list(D1 = D1, D2 = D2, D = D))

    }

    KLMatrix1 <- t(as.data.frame(lapply(bulk_1_list,function(x){
      apply(sc_1,2,function(y){KlDist(x,y)$D2})
    })))
    KL1Matrix <- max(KLMatrix1)-KLMatrix1
    rownames(KL1Matrix) <- colnames(bulk_TFscore)
    colnames(KL1Matrix) <- colnames(sc_TFscore)

    return(KL1Matrix)
  }
  simmtx <- getSimilartyMatrix(is.matrix(bulk_dataset), is.matrix(sc_dataset))
  cat("Similarty Matrix Done","\n")

  ####cell Pre-filtered
  ##X is all all single cell similarty matrix
  getCandidateCell <- function(simmtx, phenotype, screenRatio = 0.8,
                               family = "binomial"){
    X <- simmtx
    Y <- phenotype
    ratio <- screenRatio
    #    methods <- methods
    filter_num <-  dim(X)[2]*ratio
    if(family == "binomial"){
      fit <- sbisis(X = X, Y = Y,
                    candidate = filter_num,
                    method = "SBI-SIS-Pvalue",R=10)
    }
    if(family == "gaussian"){
      fit <- sbisis(X = X, Y = Y,
                    candidate = filter_num,
                    method = "SBI-SIS-Pvalue",R=10)
    }
    if(family == "cox"){
      fit <- sbisis.surv(X = X, Y = Y,
                         candidate = filter_num)
    }
    candidate_cell <- fit$candidate.set

    return(candidate_cell)
  }
  cand_cell <- getCandidateCell(simmtx, phenotype, screenRatio=screenRatio,
                                family = family)

  cat("Pre-screening Done","\n")

  ##X is all all single cell similarty matrix,
  ##Y: phenotype label,
  ##counts: the matrix of all cell
  ##sc_peak: the peaks of single cell data
  ##candidatecell: the candidate cells by runing the SBISIS
  ##ratio: The highest ratio of cells were identified
  ##resolution: the parameter findclusters()
  identifyCellSub <- function (X, Y, family = c("gaussian", "binomial", "cox")[2], object,
                               candidatecell, cutoff = 0.1, resolution = 1, group = NULL)
  {
    ratio <- cutoff
    seurat_object <- object
    seurat_object <- CreateSeuratObject(counts = as.matrix(seurat_object), assay = "RNA",
                                      min.cells = 1)
    seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
    seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = F)
    seurat_object <- ScaleData(seurat_object,verbose = F)
    seurat_object <- RunPCA(seurat_object, verbose = F)
    seurat_object <- FindNeighbors(seurat_object, dims = 1:30, verbose = F)
    seurat_object <- FindClusters(seurat_object,
                              resolution = 1)
    seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = F)

    if(!is.null(group)){

      index <- group[candidatecell]###

    } else {

      index <- seurat_object$seurat_clusters[candidatecell]
    }

    inpute_x <- X[, candidatecell]

    if(family == "cox"){
      data = list(x = inpute_x, time = Y[,1], status = Y[,2])
      fit <- SGL::SGL(data, index, type = "cox", min.frac = 0.001,
                      nlam = 1000, standardize = F, alpha = 0, lambdas = seq(0.001,1,0.001))

    }
    if(family == "gaussian"){
      data = list(x= inpute_x, y = Y)
      fit <- SGL::SGL(data, index, type = "linear", min.frac = 0.001,
                      nlam = 1000, standardize = F, alpha = 0, lambdas = seq(0.001,1,0.001))

    }
    if(family == "binomial"){
      data = list(x= inpute_x, y = Y)
      fit <- SGL::SGL(data, index, type = "logit", min.frac = 0.001,
                      nlam = 1000, standardize = F, alpha = 0, lambdas = seq(0.001,1,0.001))

    }

    for (k in 1:dim(fit$beta)[2]) {
      nonzero_ratio <- (length(which((fit$beta[, k] > 0))))/(dim(X)[2])
      num <- k
      if (nonzero_ratio < ratio)
        break
    }
    sgl_results <- fit$beta[, num]
    label_sgl <- rep("Background", dim(X)[2])
    names(label_sgl) <- colnames(X)
    label_sgl[colnames(inpute_x)[which(sgl_results > 0)]] = "PACells"
    cat("Label Summary (B/P): ", table(label_sgl))
    sc_data <- seurat_object
    sc_data@meta.data$PACells_label <- factor(label_sgl, levels = c("PACells",
                                                                    "Background"))
    return(sc_data)
  }

  sc_data <-  identifyCellSub(X=simmtx, Y=phenotype, family = c("gaussian", "binomial", "cox")[2],
                              object = sc_dataset, candidatecell=cand_cell,
                              cutoff = cutoff, resolution = 1)

  return(sc_data)
}


##  Pre-screening candidate cells based on Bcor.
##
##' @param Y a numeric matirx(phenotype).
##' @param X a numeric matirx(similarity matrix).
##' @param candidate size of candidate set.
##' @param method Method for screening procedure.
##' @param R permutation Time.
##' @import stats
##' @import utils
##' @import snowfall
##' @import gam
##' @export
sbisis <- function(Y, X, candidate = c("large"),
                   method = "SBI-SIS-Pvalue",
                   parms = list(d1=5, d2=5, df=3), R = 199, parallel = FALSE, ncore = NULL)
{
  if(parallel) {
    library(snowfall)
    sfInit(parallel = TRUE, cpus = ncore)
    sfLibrary(PACells)
  }
  n <- dim(X)[1]
  p <- dim(X)[2]
  ids <- 1:p
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  colnames(X) = paste0("X", 1:p)
  colnames(Y) = paste0("Y", 1:ncol(Y))
  if(any(apply(Y,2,anyNA))) {stop("NA appear in matrix Y")}
  if(any(apply(X,2,anyNA))) {stop("NA appear in matrix X")}

  # decide candicate size
  d_logn <- round(n/log(n))
  d <- n
  d1 <- parms$d1
  d2 <- parms$d2
  df <- parms$df
  if(is.numeric(candidate)) {
    final_d <- candidate
  } else {
    if(candidate == "small"){
      final_d <- d_logn
    } else {
      final_d <- d
    }
  }

  rcory_result <- c()

  if(method != "Interaction") {
    if(method == "SBI-SIS-Pvalue") {
      rcory_result = lapply(as.list(1:length(ids)),function(id, x = X, y = Y){
        xo <- x[,id]
        BI(y, xo, R = R)
      })
      rcory_result <- 1 - unlist(rcory_result)
    } else {
      if(ncol(Y)>1) {
        if(parallel) {
          sfExport(list = c("X", "Y"), local = TRUE)
          rcory_result <- sfLapply(as.list(1:length(ids)),function(id,x=X,y=Y){
            xo=x[,id]
            SBI(y, xo)
          })
        } else {
          rcory_result <- lapply(as.list(1:length(ids)),function(id,x=X,y=Y){
            xo=x[,id]
            SBI(y,xo)
          })
        }
      } else {
        rcory_result <- lapply(as.list(1:length(ids)),function(id, x = X, y = Y){
          xo=x[,id]
          SBI(y, xo, fast = TRUE)
        })
      }
      rcory_result <- unlist(rcory_result)
    }
    max_ids <- order(rcory_result, decreasing = T)

    method_sub <- paste0(head(unlist(strsplit(method,"-")),2),collapse='-')
    method_sub1 <- tail(unlist(strsplit(method,"-")),1)
    if(method_sub == "SBI-IISIS"){
      chooseids <- max_ids[1:d1]
      Xhavepickout <- ids[chooseids]
      Xlastpickout <- ids[chooseids]  # flag the pick out variables, used to remove the effect of selected variables

      ids <- ids[-chooseids]
      #iteration round
      if(method_sub1 == 'lm'){
        while(length(Xhavepickout) < final_d)
        {
          # lm fit for X
          Xnew <- lm(X[,ids]~X[,Xhavepickout])$resid
          # lm fit for Y
          Y <- lm(Y~X[,Xlastpickout])$resid

          # SBI-screening
          if(parallel) {
            sfExport(list = c("Xnew", "Y"))
            rcory_result <- sfLapply(as.list(1:length(ids)),function(id,x=Xnew,y=Y){
              xo=x[,id]
              SBI(y, xo)
            })
          } else {
            rcory_result <- lapply(as.list(1:length(ids)),function(id,x=Xnew,y=Y){
              xo=x[,id]
              SBI(y,xo)
            })
          }
          rcory_result <- unlist(rcory_result)
          max_ids <- order(rcory_result,decreasing=T)

          # prepare for next iteration
          chooseids <- max_ids[1:d2]
          Xhavepickout <- c(Xhavepickout,ids[chooseids])
          Xlastpickout <- ids[chooseids]
          ids <- ids[-chooseids]
        }
      }
      if(method_sub1=='gam'){
        while(length(Xhavepickout) < final_d)
        {
          # gam fit for X
          lastpickout_formula=paste0(' + s(',colnames(X)[Xlastpickout],collapse = paste0(",df = ",df,")"))
          lastpickout_formula=paste0(lastpickout_formula,paste0(",df = ",df,")"),collapse = "")
          lastpickout_dat=X[,Xlastpickout]
          Xnew=sapply(ids,function(x){
            formula_one=paste0(colnames(X)[x],"~",lastpickout_formula)
            formula_one=as.formula(formula_one)
            dat=as.data.frame(cbind(X[,x],lastpickout_dat))
            colnames(dat)[1]=colnames(X)[x]
            # dat=as.data.frame(dat)
            # colnames(dat)=paste0("X",c(x,Xhavepickout))
            gam(formula_one,data = dat)$residuals
          })

          # gam fit for Y
          dat=data.frame(Y,lastpickout_dat)
          names(dat)[1] = c("Y")
          formula_Y <- as.formula(paste("Y~", lastpickout_formula))
          Y <- gam(formula = formula_Y,data = dat)$residuals

          # SBI-screening
          rcory_result <- lapply(as.list(1:length(ids)),function(id, x = Xnew, y = Y){
            xo=x[,id]
            SBI(y, xo)
          })
          rcory_result=unlist(rcory_result)
          max_ids=order(rcory_result,decreasing=T)

          # prepare for next iteration
          chooseids <- max_ids[1:d2] #
          Xhavepickout <- c(Xhavepickout,ids[chooseids])
          Xlastpickout <- ids[chooseids]
          ids <- ids[-chooseids]
        }
      }
    } else {
      chooseids <- max_ids[1:final_d]
      Xhavepickout <- ids[chooseids]
    }
  } else {
    bcorValue <- apply(X, 2, SBI, y = Y)
    bcor2Value <- apply((X)^2, 2, SBI, y = Y)
    max_ids <- order(bcorValue, decreasing = TRUE)
    chooseids <- max_ids[1:final_d]
    Xhavepickout <- ids[chooseids]
    max_ids <- order(bcor2Value, decreasing = TRUE)
    chooseids <- max_ids[1:final_d]
    Xhavepickout <- unique(c(ids[chooseids], Xhavepickout))
  }

  return(list(candidate.set = Xhavepickout,
              candidate.data = list(X=X[,Xhavepickout],Y=Y),
              candidate.size = length(Xhavepickout)))
}


#' @param x a numeric matirx or vector
#' @param y a numeric matirx or vector
#' @param R permutation Time
#' @param seed random seed
#' @useDynLib PACells
BI <- function(x, y, R = 0, seed = 2015){
  weight <- FALSE
  x <- as.matrix(x)
  y <- as.matrix(y)
  dim_x <- dim(x)
  dim_y <- dim(y)
  n <- dim_x[1]
  p <- dim_x[2]
  q <- dim(y)[2]
  RCT <- numeric(1)
  Dx <- numeric(n*n)
  Dy <- numeric(n*n)
  dst <- TRUE
  Dx <- .C("distance", as.double(t(x)), as.double(t(Dx)), as.integer(n),as.integer(p))
  x <- matrix(Dx[[2]],n,n)
  Dy <- .C("distance", as.double(t(y)), as.double(t(Dy)), as.integer(n),as.integer(q))
  y <- matrix(Dy[[2]],n,n)
  RCT<-.C("BI", as.double(t(x)), as.double(t(y)), as.integer(n), as.integer(p), as.integer(q), as.integer(dst), HRC=as.double(RCT), as.integer(seed), as.integer(R), as.integer(weight))
  return(RCT$HRC)
}


#' @param x a numeric vector
#' @param y a numeric vector
#' @param R permutation Time
#' @param seed random seed
#' @useDynLib PACells
UBI <- function(x, y, R = 0, seed = 2015){
  weight=FALSE
  if(is.vector(x) & is.vector(y)) {
    n <- length(x)
    RCT <- numeric(1)
    RCT <- .C("UBI", as.double(x), as.double(y), as.integer(n), HRC=as.double(RCT), as.integer(seed), as.integer(R), as.integer(weight))
    return(RCT$HRC)
  }
  else {
    stop("x and y must be vector!")
  }
}


#  Calculate Ball Information statistic
#
#' @aliases SBI
#' @param x a numeric matirx
#' @param y a numeric matirx
#' @param fast Fast version for univariate case
#' @export
#' @examples
#' n <- 100
#' x <- rnorm(n)
#' y <- rnorm(n)
#' SBI(x, y)
SBI <- function(x, y, fast = FALSE){
  if(fast) {
    return(sqrt(UBI(y, x, R = 0)/sqrt(UBI(x, x, R = 0)*UBI(y, y, R = 0))))
  } else {
    return(sqrt(BI(y, x, R = 0)/sqrt(BI(x, x, R = 0)*BI(y, y, R = 0))))
  }
}


##' @param Y a numeric matirx(first column should be event time, second column should be survival status) or Surv object
##' @param X a numeric matirx.
##' @param candidate size of candidate set.
##' @param standized allows the user to standardize the covariate
##' @import stats
##' @import utils
##' @import survival
##' @export
sbisis.surv <- function(Y, X, candidate=c("large"), standized = TRUE){

  n=dim(X)[1]; p=dim(X)[2]
  n = as.numeric(n); p = as.numeric(p)
  ids=1:p
  Y=as.matrix(Y)
  X=as.matrix(X)
  colnames(X)=paste0("X",1:p)
  colnames(Y)=paste0("Y",1:ncol(Y))
  if(any(apply(Y,2,anyNA))) {stop("NA appear in matrix Y")}
  if(any(apply(X,2,anyNA))) {stop("NA appear in matrix X")}

  # decide candicate size
  d_logn=round(n/log(n))
  d=n
  if(is.numeric(candidate)) {
    final_d <- candidate
  } else if(candidate=="small"){
    final_d=d_logn
  } else {
    final_d=d
  }

  # prepare for screening
  time=Y[,1]
  delta=Y[,2]
  ord.t = sort(time)
  ix = order(time)
  ord.delta = delta[ix]
  xo=X[ix,]
  if(standized) {
    xo = apply(xo, 2, scale)
  }

  # SBI Screening(survival)
  fitc = survfit(Surv(time,1-delta)~1)
  Sc = fitc$surv
  if(length(unique(ord.t)) != n) {
    rep_num = as.data.frame(table(ord.t))[, "Freq"]
    Sc = mapply(function(x, y) {
      rep(x, y)
    }, Sc, rep_num, SIMPLIFY = FALSE)
    Sc= unlist(Sc)
  }
  # num_non_zero = n - sum()
  rcory_result=apply(xo,2,function(x){
    SRCT(x = x,t = ord.t,delta = ord.delta,Sc = Sc,n = n)
  })
  rcory_result=unlist(rcory_result)
  max_ids=order(rcory_result,decreasing=T)
  chooseids=max_ids[1:final_d]
  Xhavepickout=ids[chooseids]

  return(list(candidate.set=Xhavepickout,
              candidate.data=list(X=X[,Xhavepickout],Y=Y),
              candidate.size=length(Xhavepickout)))
}


# @param x ordered covariate
# @param t ordered survival event time
# @param delta ordered survival event status
# @param Sc Survfit object
# @param n Sample size
# @useDynLib PACells
# @export
#
SRCT <- function(x, t, delta, Sc, n){
  RCT <- numeric(1)
  RCT<-.C("SRCT", as.double(t(x)), as.double(t(t)), as.double(t(delta)),
          as.double(t(Sc)), as.integer(n), RC=as.double(RCT))
  return(RCT$RC)
}


# @param time Time
# @param status Status
# @param x X
# @return Survival Ball Correlation
# @import stats
# @import utils
# @import survival
# @export
# @examples
# data(survdat)
# bcor.surv(time = survdat[["time"]], status = survdat[["status"]], x = survdat[,3])
bcor.surv <- function(time, status, x) {
  n = length(time)
  ord.t = sort(time)
  ix = order(time)
  ord.delta = status[ix]
  x = scale(x[ix])
  #
  fitc = survfit(Surv(time,1-status)~1)
  Sc = fitc$surv
  if(length(unique(ord.t)) != n) {
    rep_num = as.data.frame(table(ord.t))[, "Freq"]
    Sc = mapply(function(x, y) {
      rep(x, y)
    }, Sc, rep_num, SIMPLIFY = FALSE)
    Sc= unlist(Sc)
  }
  #
  SRCT(x = x,t = ord.t, delta = ord.delta, Sc = Sc, n = n)
}

