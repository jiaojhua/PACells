#' Create a Seurat Object for chromatin accessibility data
#'
#' This function builds a Seurat object from ATAC-seq count data and is designed
#' to work for both bulk ATAC-seq and single-cell ATAC-seq.
#' Input takes a count matrix. The expected format of the input matrix is
#' features x cells. A set of genomic ranges (peaks or bins) must be supplied
#' along with the matrix, with the length of the ranges equal to the number of
#' rows in the matrix.
#'
#'Conceptually, the function performs two steps:
#' (1) create a Signac \code{ChromatinAssay} that stores the peak/bin count matrix
#' together with genomic coordinates and genome annotation; and
#' (2) wrap the assay into a Seurat object and attach sample/cell-level metadata.
#'
#' @param matrix A count matrix of chromatin accessibility in features × cells/samples format.
#' Rows correspond to genomic features (peaks or bins) and columns correspond to cells (scATAC-seq)
#' or bulk samples (bulk ATAC-seq). Sparse matrices (e.g., \code{dgCMatrix}) are recommended.
#' @param ranges A \code{GRanges} object giving the genomic coordinates of the features (peaks/bins).
#' Its length must equal the number of rows in \code{matrix}, and the order must match the row order
#' of \code{matrix} exactly.
#' @param meta.data A \code{data.frame} containing per-cell/per-sample metadata (e.g., sample ID,
#' condition, patient ID). Row names should match \code{colnames(matrix)},
#' so that each column in \code{matrix} can be correctly annotated in the Seurat object.
#' @param genome UCSC genome identifier used to label the assay genome (e.g., \code{"hg19"},
#' \code{"hg38"}). This should be consistent with the coordinate system of \code{ranges}.
#' @param min.cells Feature filtering threshold: keep only peaks/bins that are detected in at least
#' \code{min.cells} columns (cells/samples). This helps remove extremely rare features and reduce noise.
#' Default is \code{1}, meaning a feature is retained if it appears in more than one cell/sample.
#' Users can adjust this parameter according to dataset size and sparsity (e.g., using a larger value to apply stricter filtering).
#'
#' @return A Seurat object with an \code{"ATAC"} assay containing a Signac \code{ChromatinAssay}.
#' The assay stores the input count matrix and genomic ranges, and \code{meta.data} is attached as
#' cell/sample-level metadata.
#'
#' @export

creatSeurat <- function(matrix, ranges, meta.data, genome, min.cells = 1){

  object <- CreateSeuratObject(
    counts = CreateChromatinAssay(
      counts = matrix,
      ranges = ranges,
       genome = genome,
      min.cells = min.cells,
      min.features = 0),
    assay = "ATAC",
    min.cells = min.cells,
    min.features = 0,
    meta.data = meta.data
  )
  return(object)

}

#' Get TF motifs used by PACells
#'
#' This function prepares a transcription factor (TF) motif collection for PACells.
#' In PACells, TF motifs provide sequence-based regulatory patterns that are used to summarize
#' chromatin accessibility into motif-level signals, which can help characterize regulatory programs
#' and support bulk–single-cell matching in the PACells workflow.
#'
#' @param species Species identifier for querying JASPAR motifs. It can be provided as a Latin name
#' (e.g., \code{"Homo sapiens"}) or an NCBI taxonomy ID (e.g., \code{9606}). Note that, in the
#' current implementation, \code{species} is internally set to \code{"Homo sapiens"} when
#' \code{database = "JASPAR"}.
#'
#' @param database Motif source to use. Supported options are \code{"JASPAR"} and \code{"cisBP"}.
#' Default is \code{"JASPAR"}.
#'
#' @return A PFMatrixList object.
#'
#' @export

getMotifs <- function(species = "Homo sapiens",
                      database = c("JASPAR", "cisBP")[1], JASPAR = JASPAR2022){


  if(database == "JASPAR"){
    species <- species
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
#'  \code{PACells} links bulk-level clinical phenotypes (e.g., case/control status, quantitative traits,
#' prognosis, treatment response, or survival outcomes) to single-cell ATAC-seq profiles and identifies
#' a subset of cells whose chromatin accessibility patterns are most associated with the given phenotype.
#' The function returns a Seurat object in which each cell is labeled as \code{"PACells"} or
#' \code{"Background"} in the metadata.
#'
#' @param sc_dataset A Seurat object containing single-cell ATAC-seq data. It should include an ATAC assay
#' (default \code{"ATAC"}) with counts and genomic ranges.
#' @param bulk_dataset A Seurat object containing bulk ATAC-seq samples in the same feature space.
#' @param phenotype Phenotype annotation for bulk samples.
#' \itemize{
#'   \item For \code{family = "binomial"} or \code{"gaussian"}: a vector with length equal to the number of
#'   bulk samples.
#'   \item For \code{family = "cox"}: a two-column matrix/data.frame where the first column is survival time
#'   and the second column is event status, with rows aligned to bulk samples.
#' }
#' @param motifs TF motif definitions used by PACells (e.g., obtained from \code{getMotifs()}), provided as
#' a motif collection compatible with the motif processing steps in the pipeline.
#' @param cutoff Cutoff for the percentage of the PACells selected cells in total cells. This parameter is used to
#' restrict the number of the PACells selected cells. A smaller cutoff value (default \code{10\%}) is recommended
#' depending on the input data. Users can adjust it based on dataset size and the desired specificity.
#' @param screenRatio Pre-screen candidate cells that are strongly related to the phenotype.
#' This parameter is used to pre-screen the percentage of candidate cells based Bcor measure.
#' If a numeric value is provided, it directly indicates the
#' proportion of cells to keep (default \code{0.80}). If \code{"auto"} is used, PACells will select
#' a suitable proportion by checking whether the pre-screening results remain stable under repeated
#' subsampling of bulk samples (controlled by auto_subsample and auto_stab_threshold).
#' @param family Response type for the regression model. It depends on the type of the given phenotype.
#' \code{"gaussian"} for continuous outcomes, \code{"binomial"} for binary outcomes, and \code{"cox"} for
#' survival outcomes.
#' @param method Method for calculating the similarity matrix of bulk and single cells.
#' Options include \code{"KL"}, \code{"Pearson"}, and \code{"Spearman"}. The default is KL divergence.
#' @param level Logical value. If TRUE, calculating similarity based on an overall activity level profile.
#' If FALSE, calculating similarity based on relative activity patterns between features. The default is TRUE.
#' @param sc_refgenome Reference genome for sc_dataset ("hg38" or "hg19"), used to
#' interpret genomic coordinates during motif-related processing.
#' @param bulk_refgenome A reference genome for bulk ATAC-seq data("hg38" or "hg19"), used to
#' interpret genomic coordinates during motif-related processing.
#' @param res the resolution parameter (default 1) in the FindClusters function, which is used to group cells.
#' @param dims_Neighbors Dimensions of reduction to construct nearest-neighbor graph.
#' @param dims_UMAP Dimensions used to compute UMAP embedding.
#' @param group A vector of grouping information of single cells, when is provided by the user, is not calculated by the Signac package.
#' The vector should be aligned with cells in \code{sc_dataset} (by order or by names).
#' @param auto_subsample Subsampling proportion of bulk samples used by the adaptive
#'   candidate pre-screening strategy when \code{screenRatio = "auto"}.
#'   The default is \code{0.7}.
#' @param auto_stab_threshold Stability threshold used by \code{screenRatio = "auto"} to select
#'   the smallest candidate ratio with stable pre-screening results across subsamples.
#'   The default is \code{0.7}.
#'   @param batch Batch correction option applied to the scATAC-seq data.
#' Supported options:
#' \itemize{
#'   \item \code{"none"}: no batch correction (default).
#'   \item \code{"harmony"}: Harmony integration on the LSI reduction.
#'   \item \code{"seurat"}: Seurat reciprocal-LSI anchor-based integration (\code{integrated_lsi}).
#' }
#' In our practice, we do not enforce batch correction by default, because many scATAC-seq datasets are
#' generated within a single study/batch and unnecessary correction may risk over-correcting and weakening
#' genuine biological signals. When batch effects are evident (e.g., UMAP shows batch-driven separation or
#' the data combine multiple studies/platforms/labs), we recommend using \code{batch = "harmony"}.
#' Across our additional comparisons, Harmony typically improved or matched identification performance
#' under pronounced cross-study heterogeneity, while keeping results broadly consistent when batch effects
#' were weak, providing a favorable balance between mitigating batch-driven structure and preserving
#' biological signals.
#' @param Batch_variable Metadata column name in \code{sc_dataset@meta.data} indicating batch/sample IDs
#' used by \code{batch = "harmony"} or \code{batch = "seurat"}. Default is \code{"Sample"}.
#' Users should ensure this column exists and correctly represents batch structure.
#' @param dim_Batch Dimensions used for batch correction (Harmony/Seurat integration).
#'
#'
#' @return A Seurat object containing PACells results. In particular, \code{sc_dataset@meta.data$PACells_label}
#' is added (levels: "PACells" and "Background"), along with any
#' results generated during the run when applicable.
#'
#'
#' @export

PACells <- function(sc_dataset, bulk_dataset, phenotype, motifs,
                    cutoff = 0.1, screenRatio = 0.8,
                    family = c("binomial", "gaussian", "cox")[1],
                    method = c("KL", "Pearson", "Spearman")[1],
                    level = TRUE,
                    sc_refgenome = c("hg38", "hg19")[1],
                    bulk_refgenome = c("hg38", "hg19")[1],
                    res = 1, dims_Neighbors = 2:30, dims_UMAP = 2:30,
                    group = NULL, auto_subsample = 0.7,
                    auto_stab_threshold = 0.7,
                    batch = c("none", "harmony", "seurat")[1],
                    Batch_variable = "Sample",
                    dim_Batch = 2:30
                    )
{


  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
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
  sc_TFscore[is.na(sc_TFscore)]=0
  bulk_TFscore[is.na(bulk_TFscore)]=0
  cat("TF Activity Matrix done","\n")

  if(method == "KL"){

    getSimilartyMatrix <- function(bulk_TFscore, sc_TFscore, level = TRUE){
      dataset0 <- cbind(bulk_TFscore,sc_TFscore)
      dataset1 <- preprocessCore::normalize.quantiles(dataset0, keep.names = TRUE)
      bulk_TFscore1 <- dataset1[,1:dim(bulk_TFscore)[2]]
      sc_TFscore1 <- dataset1[,-c(1:dim(bulk_TFscore)[2])]
      if (level){
      nbins <- seq(min(dataset1), max(dataset1),length.out = dim(dataset1)[1]*2)
      sc_TFscore_count <- apply(sc_TFscore1, 2, function(x){hist(x, breaks = nbins, plot = FALSE)$counts})
      bulk_TFscore_count <- apply(bulk_TFscore1, 2, function(x){hist(x, breaks = nbins, plot = FALSE)$counts})
      } else {
        bulk_TFscore_count <- apply(bulk_TFscore1, 2, robust_minmax)
        sc_TFscore_count   <- apply(sc_TFscore1,   2, robust_minmax)
      }
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

    simmtx <- getSimilartyMatrix(bulk_TFscore, sc_TFscore, level)
  }

  if(method == "Pearson"){
    getPearsonMatrix <- function(bulk_TFscore, sc_TFscore){
      dataset0 <- cbind(bulk_TFscore,sc_TFscore)
      dataset1 <- preprocessCore::normalize.quantiles(dataset0, keep.names = TRUE)
      bulk_TFscore1 <- dataset1[,1:dim(bulk_TFscore)[2]]
      sc_TFscore1 <- dataset1[,-c(1:dim(bulk_TFscore)[2])]
      pearsonMatrix <- cor(bulk_TFscore1, sc_TFscore1, method = c("pearson", "spearman")[1])

      return(pearsonMatrix)
    }
    simmtx <- getPearsonMatrix(bulk_TFscore, sc_TFscore)
  }

  if(method == "Spearman"){
    getSpearmanMatrix <- function(bulk_TFscore, sc_TFscore){
      dataset0 <- cbind(bulk_TFscore,sc_TFscore)
      dataset1 <- preprocessCore::normalize.quantiles(dataset0, keep.names = TRUE)
      bulk_TFscore1 <- dataset1[,1:dim(bulk_TFscore)[2]]
      sc_TFscore1 <- dataset1[,-c(1:dim(bulk_TFscore)[2])]
      spearsonMatrix <- cor(bulk_TFscore1, sc_TFscore1, method = c("pearson", "spearman")[2])

      return(spearsonMatrix)
    }
    simmtx <- getSpearmanMatrix(bulk_TFscore, sc_TFscore)
  }

  cat("Similarty Matrix Done","\n")

  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  getCandidateCell <- function(simmtx, phenotype, screenRatio = 0.8,
                               family = "binomial",
                               auto_ratios = c(0.7, 0.8, 0.9),
                               auto_B = 10, auto_subsample = 0.7,
                               auto_stab_threshold = 0.7){

    X <- simmtx
    Y <- phenotype

    jaccard <- function(a, b) {
      a <- unique(a); b <- unique(b)
      length(intersect(a, b)) / length(union(a, b))
    }

    run_one_screen <- function(Xb, Yb, ratio) {
      m <- dim(Xb)[2]
      filter_num <- max(1, (m * ratio))

      if(family == "binomial"){
        fit <- sbisis(X = Xb, Y = Yb,
                      candidate = filter_num,
                      method = "SBI-SIS-Pvalue",R=10)
      }
      if(family == "gaussian"){
        fit <- sbisis(X = Xb, Y = Yb,
                      candidate = filter_num,
                      method = "SBI-SIS-Pvalue",R=10)
      }
      if(family == "cox"){
        fit <- sbisis.surv(X = Xb, Y = Yb,
                           candidate = filter_num)
      }
      candidate_cell <- fit$candidate.set
    }

    autoPickRatio <- function(X, Y) {
      n <- nrow(X)

      stab_scores <- sapply(auto_ratios, function(r) {
        sets <- vector("list", auto_B)

        for (b in seq_len(auto_B)) {
          set.seed(10000 + b)
          idx <- sample(seq_len(n), size = max(2, floor(n * auto_subsample)))

          Xb <- X[idx, , drop = FALSE]
          if (family == "cox") {
            Yb <- Y[idx, , drop = FALSE]
          } else {
            Yb <- Y[idx]
          }

          sets[[b]] <- run_one_screen(Xb, Yb, r)
        }


        jac <- c()
        if (auto_B >= 2) {
          for (i in 1:(auto_B - 1)) {
            for (j in (i + 1):auto_B) {
              jac <- c(jac, jaccard(sets[[i]], sets[[j]]))
            }
          }
        }
        mean(jac, na.rm = TRUE)
      })

      ok <- which(stab_scores >= auto_stab_threshold)
      if (length(ok) > 0) {
        min(auto_ratios[ok])
      } else {
        0.8
      }
    }

    if (is.character(screenRatio) && screenRatio == "auto") {
      ratio <- autoPickRatio(X, Y)
    } else {
      ratio <- screenRatio
    }


    filter_num <- max(1, (dim(X)[2] * ratio))

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
  }

  cand_cell <- getCandidateCell(simmtx, phenotype, screenRatio=screenRatio,
                                family = family,
                                auto_subsample = auto_subsample,
                                auto_stab_threshold = auto_stab_threshold)

  cat("Pre-screening Done","\n")

  ##X is all all single cell similarty matrix,
  ##Y: phenotype label,
  ##counts: the matrix of all cell
  ##sc_peak: the peaks of single cell data
  ##candidatecell: the candidate cells by runing the SBISIS
  ##ratio: The highest ratio of cells were identified
  ##resolution: the parameter findclusters()
  identifyCellSub <- function (X, Y, family = c("binomial", "gaussian", "cox")[1], object,
                               candidatecell, cutoff = 0.1, resolution = 1,
                               dims_Neighbors = 2:30, dims_UMAP = 2:30, group = NULL,
                               batch = c("none", "harmony", "seurat")[1],
                               Batch_variable = "Sample",
                               dim_Batch = 2:30 )
  {
    ratio <- cutoff
    seurat_object <- object
    seurat_object <- RunTFIDF(seurat_object, verbose = F)
    seurat_object <- FindTopFeatures(seurat_object, min.cutoff = "q10",
                                     verbose = F)
    seurat_object <- RunSVD(seurat_object, verbose = F)

    if (batch == "none") {
      reduction_use <- "lsi"
    }

    if(batch == "harmony") {
      seurat_object <- RunHarmony(
        object = seurat_object,
        max_iter = 10,
        group.by.vars = Batch_variable,
        dims.use = dim_Batch,
        reduction = 'lsi',
        assay.use = 'ATAC',
        project.dim = FALSE,
        verbose = F
        )
      reduction_use <- "harmony"
    }

    if (batch == "seurat") {
      obj_list <- SplitObject(seurat_object, split.by = Batch_variable)

      obj_list <- lapply(
        obj_list,
        function(x) {
          if ("ATAC" %in% names(x@assays)) {
            DefaultAssay(x) <- "ATAC"
          }
          x <- FindTopFeatures(x, min.cutoff = "q10", verbose = FALSE)
          x <- RunTFIDF(x, verbose = FALSE)
          x <- RunSVD(x, verbose = FALSE)
          return(x)
        }
      )

      anchors <- FindIntegrationAnchors(
        object.list    = obj_list,
        anchor.features = rownames(obj_list[[1]]),
        reduction      = "rlsi",
        dims           = dim_Batch
      )

      seurat_object <- IntegrateEmbeddings(
        anchorset         = anchors,
        reductions         = seurat_object[["lsi"]],
        new.reduction.name = "integrated_lsi",
        dims.to.integrate = c(1, dim_Batch)
      )

      reduction_use <- "integrated_lsi"

    }

    seurat_object <- RunUMAP(object = seurat_object, reduction = reduction_use,
                             dims = dims_UMAP, verbose = F)
    seurat_object <- FindNeighbors(object = seurat_object, reduction = reduction_use,
                                   dims = dims_Neighbors, verbose = F)
    seurat_object <- FindClusters(object = seurat_object, verbose = FALSE,
                                  algorithm = 3, resolution = resolution)
    if(!is.null(group)){

      index <- group[candidatecell]

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

  sc_data <-  identifyCellSub(X=simmtx, Y=phenotype, family = family,
                              object = sc_dataset, candidatecell=cand_cell,
                              cutoff = cutoff, resolution = res,
                              dims_Neighbors = dims_Neighbors, dims_UMAP = dims_UMAP,
                              group=group, batch = batch, Batch_variable = Batch_variable,
                              dim_Batch = dim_Batch)

  return(sc_data)
}

#' PPACells.Covar: identify clinical phenotype-associated cell states from scATAC-seq data with covariates
#'
#'  \code{PACells.Covar} extends \code{PACells} by allowing users to include bulk-sample covariates when linking
#' a bulk phenotype to single cells. This function is intended for practical scenarios where the phenotype of
#' interest may be correlated with known sample-level variables (e.g., demographic variables, clinical factors,
#' available for bulk samples). By providing \code{covariates}, users can evaluate
#' phenotype-associated cells while accounting for these additional variables.
#'
#' The function returns the input scATAC-seq Seurat object with an added metadata column
#' \code{PACells_label}, which labels each cell as \code{"PACells"} (selected phenotype-associated cells) or
#' \code{"Background"} (all other cells). In addition, the function may add intermediate dimensional reductions
#' and embedding (e.g., LSI, Harmony/integrated LSI, UMAP) used during the analysis.
#'
#'#' \strong{Input alignment is important:}
#' \itemize{
#'   \item \code{phenotype} and \code{covariates} should follow the same bulk-sample order as in \code{bulk_dataset}.
#'   \item For \code{family = "cox"}, \code{phenotype} should contain two columns (\code{time} and \code{status}).
#' }
#'
#' @param sc_dataset A Seurat object containing scATAC-seq data. It should include an \code{"ATAC"} assay with
#' peak/bin counts and genomic ranges.
#' @param bulk_dataset bulk_dataset A Seurat object containing bulk ATAC-seq data. Each column represents one bulk sample.
#' Bulk samples should correspond to entries in \code{phenotype} and rows in \code{covariates}.
#' @param phenotype Phenotype annotation for bulk samples.
#' \itemize{
#'   \item For \code{family = "binomial"} and \code{family = "gaussian"}: a vector aligned to bulk samples.
#'   \item For \code{family = "cox"}: a two-column matrix/data.frame with survival \code{time} and event
#'   \code{status}, aligned to bulk samples.
#' }
#' @param covariates A data.frame of covariates for bulk samples. Rows should correspond to bulk samples in the
#' same order as \code{phenotype}. Columns are covariate variables. Users may include clinical
#' covariates that they would like to account for in the phenotype association step.
#' @param motifs TF motif definitions (e.g., from \code{getMotifs()}) used in the PACells workflow.
#' @param cutoff Maximum proportion of cells labeled as \code{"PACells"} among all single cells.
#' This parameter restricts the final number of selected cells.
#' A smaller cutoff value (default \code{10\%}) is recommended depending on the input data.
#' Users can adjust it depending on data size and how stringent the cell selection should be.
#'
#' @param screenRatio Pre-screen candidate cells that are strongly related to the phenotype.
#' This parameter controls the percentage of cells retained as candidates before the final labeling step.
#' If a numeric value is provided, it directly specifies the proportion of cells to keep (default \code{0.80}).
#' If \code{screenRatio = "auto"}, PACells selects a suitable proportion by checking the stability of the
#' pre-screening results under repeated subsampling of bulk samples (controlled by \code{auto_subsample} and
#' \code{auto_stab_threshold}).
#' @param family Response type for the regression model. It depends on the type of the given phenotype and
#' can be \code{family = gaussian} for linear regression, \code{family = binomial} for classification,
#' or \code{family = cox} for Cox regression.
#' @param method Method for calculating the similarity matrix of bulk and single cells.
#' Options include \code{"KL"}, \code{"Pearson"}, and \code{"Spearman"}. The default is KL divergence.
#' @param level Logical value. If TRUE, calculating similarity based on an overall activity level profile.
#' If FALSE, calculating similarity based on relative activity patterns between features. The default is TRUE.
#' @param sc_refgenome Reference genome for the scATAC-seq object ("hg38" or "hg19"). This should
#' match the coordinate system of the peak ranges stored in \code{sc_dataset}.
#' @param sc_refgenome Reference genome for the bulk ATAC-seq object ("hg38" or "hg19"). This should
#' match the coordinate system of the peak ranges stored in \code{bulk_dataset}.
#' @param res the resolution parameter (default 1) in the FindClusters function, which is used to group cells.
#' @param dims_Neighbors Dimensions of reduction to construct nearest-neighbor graph.
#' @param dims_UMAP Dimensions of reduction to UMAP.
#' @param group Optional user-provided grouping vector for single cells. A vector of grouping information of
#' single cells, when is provided by the user, this grouping is used instead of clusters computed internally.
#'@param auto_subsample Subsampling proportion of bulk samples used when \code{screenRatio = "auto"}.
#' Default is \code{0.7}.#'
#' @param auto_stab_threshold Stability threshold used when \code{screenRatio = "auto"}.
#' Default is \code{0.7}.#'
#' @param batch Batch correction option applied to the scATAC-seq preprocessing step.
#' Supported options:
#' \itemize{
#'   \item \code{"none"}: no batch correction (default).
#'   \item \code{"harmony"}: Harmony integration on the LSI reduction.
#'   \item \code{"seurat"}: Seurat reciprocal-LSI anchor-based integration.
#' }
#' In our practice, we do not enforce batch correction by default, because many scATAC-seq datasets are generated
#' within a single study and unnecessary correction may risk over-correcting. When batch effects are evident
#' (e.g., cells separate by batch/study in uncorrected embeddings or data are combined across studies/platforms),
#' we recommend using \code{batch = "harmony"}, which in our additional comparisons generally improved reliability
#' under pronounced cross-study heterogeneity while remaining consistent when batch effects were weak.
#' @param Batch_variable Metadata column name in \code{sc_dataset@meta.data} indicating batch/sample IDs used by
#' \code{batch = "harmony"} or \code{batch = "seurat"}. Default is \code{"Sample"}.
#' @param dim_Batch Dimensions used for batch correction (Harmony/Seurat integration).
#'
#' @return This function returns a Seurat-class object. It contains the results that cell state identified by
#' PACells as strongly associated with the clinical phenotype.
#' A Seurat object based on \code{sc_dataset} with an additional metadata column PACells_label
#' (levels: "PACells" and "Background"). Cells labeled "PACells" are those identified as most
#' associated with the bulk phenotype after accounting for the provided covariates.
#'
#' @export

PACells.Covar <- function(sc_dataset, bulk_dataset, phenotype,covariates, motifs,
                          cutoff = 0.1, screenRatio = 0.8,
                          family = c("binomial", "gaussian", "cox")[1],
                          method = c("KL", "Pearson", "Spearman")[1],
                          level = TRUE,
                          sc_refgenome = c("hg38", "hg19")[1],
                          bulk_refgenome = c("hg38", "hg19")[1],
                          res = 1, dims_Neighbors = 2:30, dims_UMAP = 2:30,
                          group = NULL, auto_subsample = 0.7,
                          auto_stab_threshold = 0.7,
                          batch = c("none", "harmony", "seurat")[1],
                          Batch_variable = "Sample",
                          dim_Batch = 2:30)
{

  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
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
  sc_TFscore[is.na(sc_TFscore)]=0
  bulk_TFscore[is.na(bulk_TFscore)]=0
  cat("TF Activity Matrix done","\n")

    if(method == "KL"){

      getSimilartyMatrix <- function(bulk_TFscore, sc_TFscore, level = TRUE){
        dataset0 <- cbind(bulk_TFscore,sc_TFscore)
        dataset1 <- preprocessCore::normalize.quantiles(dataset0, keep.names = TRUE)
        bulk_TFscore1 <- dataset1[,1:dim(bulk_TFscore)[2]]
        sc_TFscore1 <- dataset1[,-c(1:dim(bulk_TFscore)[2])]
        if (level){
          nbins <- seq(min(dataset1), max(dataset1),length.out = dim(dataset1)[1]*2)
          sc_TFscore_count <- apply(sc_TFscore1, 2, function(x){hist(x, breaks = nbins, plot = FALSE)$counts})
          bulk_TFscore_count <- apply(bulk_TFscore1, 2, function(x){hist(x, breaks = nbins, plot = FALSE)$counts})
        } else {
          bulk_TFscore_count <- apply(bulk_TFscore1, 2, robust_minmax)
          sc_TFscore_count   <- apply(sc_TFscore1,   2, robust_minmax)
        }
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

      simmtx <- getSimilartyMatrix(bulk_TFscore, sc_TFscore, level)
    }


  if(method == "Pearson"){
    getPearsonMatrix <- function(bulk_TFscore, sc_TFscore){
      dataset0 <- cbind(bulk_TFscore,sc_TFscore)
      dataset1 <- preprocessCore::normalize.quantiles(dataset0, keep.names = TRUE)
      bulk_TFscore1 <- dataset1[,1:dim(bulk_TFscore)[2]]
      sc_TFscore1 <- dataset1[,-c(1:dim(bulk_TFscore)[2])]
      pearsonMatrix <- cor(bulk_TFscore1, sc_TFscore1, method = c("pearson", "spearman")[1])

      return(pearsonMatrix)
    }

    simmtx <- getPearsonMatrix(bulk_TFscore, sc_TFscore)

  }


  if(method == "Spearman"){
    getSpearmanMatrix <- function(bulk_TFscore, sc_TFscore){
      dataset0 <- cbind(bulk_TFscore,sc_TFscore)
      dataset1 <- preprocessCore::normalize.quantiles(dataset0, keep.names = TRUE)
      bulk_TFscore1 <- dataset1[,1:dim(bulk_TFscore)[2]]
      sc_TFscore1 <- dataset1[,-c(1:dim(bulk_TFscore)[2])]
      spearsonMatrix <- cor(bulk_TFscore1, sc_TFscore1, method = c("pearson", "spearman")[2])

      return(spearsonMatrix)
    }

    simmtx <- getSpearmanMatrix(bulk_TFscore, sc_TFscore)

  }

  cat("Similarty Matrix Done","\n")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  getCandidateCell <- function(simmtx, phenotype, screenRatio = 0.8,
                               family = "binomial",
                               auto_ratios = c(0.7, 0.8, 0.9),
                               auto_B = 10, auto_subsample = 0.7,
                               auto_stab_threshold = 0.7){

    X <- simmtx
    Y <- phenotype

    jaccard <- function(a, b) {
      a <- unique(a); b <- unique(b)
      length(intersect(a, b)) / length(union(a, b))
    }

    run_one_screen <- function(Xb, Yb, ratio) {
      m <- dim(Xb)[2]
      filter_num <- max(1, (m * ratio))

      if(family == "binomial"){
        fit <- sbisis(X = Xb, Y = Yb,
                      candidate = filter_num,
                      method = "SBI-SIS-Pvalue",R=10)
      }
      if(family == "gaussian"){
        fit <- sbisis(X = Xb, Y = Yb,
                      candidate = filter_num,
                      method = "SBI-SIS-Pvalue",R=10)
      }
      if(family == "cox"){
        fit <- sbisis.surv(X = Xb, Y = Yb,
                           candidate = filter_num)
      }
      candidate_cell <- fit$candidate.set
    }

    autoPickRatio <- function(X, Y) {
      n <- nrow(X)

      stab_scores <- sapply(auto_ratios, function(r) {
        sets <- vector("list", auto_B)

        for (b in seq_len(auto_B)) {
          set.seed(10000 + b)
          idx <- sample(seq_len(n), size = max(2, floor(n * auto_subsample)))

          Xb <- X[idx, , drop = FALSE]
          if (family == "cox") {
            Yb <- Y[idx, , drop = FALSE]
          } else {
            Yb <- Y[idx]
          }

          sets[[b]] <- run_one_screen(Xb, Yb, r)
        }


        jac <- c()
        if (auto_B >= 2) {
          for (i in 1:(auto_B - 1)) {
            for (j in (i + 1):auto_B) {
              jac <- c(jac, jaccard(sets[[i]], sets[[j]]))
            }
          }
        }
        mean(jac, na.rm = TRUE)
      })

      ok <- which(stab_scores >= auto_stab_threshold)
      if (length(ok) > 0) {
        min(auto_ratios[ok])
      } else {
        0.8
      }
    }

    if (is.character(screenRatio) && screenRatio == "auto") {
      ratio <- autoPickRatio(X, Y)
    } else {
      ratio <- screenRatio
    }


    filter_num <- max(1, (dim(X)[2] * ratio))

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
  }
  cand_cell <- getCandidateCell(simmtx, phenotype, screenRatio=screenRatio,
                                family = family,
                                auto_subsample = auto_subsample,
                                auto_stab_threshold = auto_stab_threshold)

  cat("Pre-screening Done","\n")

  ##X is all all single cell similarty matrix,
  ##Y: phenotype label,
  ##counts: the matrix of all cell
  ##sc_peak: the peaks of single cell data
  ##candidatecell: the candidate cells by runing the SBISIS
  ##ratio: The highest ratio of cells were identified
  ##resolution: the parameter findclusters()
  identifyCellSub.Covar <- function (X, Y, covariates,family = c("binomial", "gaussian", "cox")[1], object,
                                     candidatecell, cutoff = 0.1, resolution = 1,
                                     dims_Neighbors = 2:30, dims_UMAP = 2:30, group = NULL,
                                     batch = c("none", "harmony", "seurat")[1],
                                     Batch_variable = "Sample",
                                     dim_Batch = 2:30)
  {
    ratio <- cutoff
    seurat_object <- object
    seurat_object <- RunTFIDF(seurat_object, verbose = F)
    seurat_object <- FindTopFeatures(seurat_object, min.cutoff = "q10",
                                     verbose = F)
    seurat_object <- RunSVD(seurat_object, verbose = F)


    if (batch == "none") {
      reduction_use <- "lsi"
    }

    if(batch == "harmony") {
      seurat_object <- RunHarmony(
        object = seurat_object,
        max_iter = 10,
        group.by.vars = Batch_variable,
        dims.use = dim_Batch,
        reduction = 'lsi',
        assay.use = 'ATAC',
        project.dim = FALSE,
        verbose = F
      )
      reduction_use <- "harmony"
    }

    if (batch == "seurat") {
      obj_list <- SplitObject(seurat_object, split.by = Batch_variable)

      obj_list <- lapply(
        obj_list,
        function(x) {
          if ("ATAC" %in% names(x@assays)) {
            DefaultAssay(x) <- "ATAC"
          }
          x <- FindTopFeatures(x, min.cutoff = "q10", verbose = FALSE)
          x <- RunTFIDF(x, verbose = FALSE)
          x <- RunSVD(x, verbose = FALSE)
          return(x)
        }
      )

      anchors <- FindIntegrationAnchors(
        object.list    = obj_list,
        anchor.features = rownames(obj_list[[1]]),
        reduction      = "rlsi",
        dims           = dim_Batch
      )

      seurat_object <- IntegrateEmbeddings(
        anchorset         = anchors,
        reductions         = seurat_object[["lsi"]],
        new.reduction.name = "integrated_lsi",
        dims.to.integrate = c(1, dim_Batch)
      )

      reduction_use <- "integrated_lsi"

    }


    seurat_object <- RunUMAP(object = seurat_object, reduction = reduction_use,
                             dims = dims_UMAP, verbose = F)
    seurat_object <- FindNeighbors(object = seurat_object, reduction = reduction_use,
                                   dims = dims_Neighbors, verbose = F)
    seurat_object <- FindClusters(object = seurat_object, verbose = FALSE,
                                  algorithm = 3, resolution = resolution)
    if(!is.null(group)){

      index <- group[candidatecell]

    } else {
      index <- seurat_object$seurat_clusters[candidatecell]
    }
    index_covar <- c(rep(length(table(index))+1, dim(covariates)[2]),index)
    cell_x <- X[, candidatecell]
    inpute_x <- cbind(covariates,cell_x)

    if(family == "cox"){
      data = list(x = inpute_x, time = Y[,1], status = Y[,2])
      fit <- SGL::SGL(data, index_covar, type = "cox", min.frac = 0.001,
                      nlam = 1000, standardize = F, alpha = 0, lambdas = seq(0.001,1,0.001))

    }
    if(family == "gaussian"){
      data = list(x= inpute_x, y = Y)
      fit <- SGL::SGL(data, index_covar, type = "linear", min.frac = 0.001,
                      nlam = 1000, standardize = F, alpha = 0, lambdas = seq(0.001,1,0.001))

    }
    if(family == "binomial"){
      data = list(x= inpute_x, y = Y)
      fit <- SGL::SGL(data, index_covar, type = "logit", min.frac = 0.001,
                      nlam = 1000, standardize = F, alpha = 0, lambdas = seq(0.001,1,0.001))

    }
    beta_cell <- fit$beta[-dim(covariates)[2],]
    for (k in 1:dim(beta_cell)[2]) {
      nonzero_ratio <- (length(which((beta_cell[, k] > 0))))/(dim(X)[2])
      num <- k
      if (nonzero_ratio < ratio)
        break
    }

    sgl_results <- beta_cell[, num]
    label_sgl <- rep("Background", dim(X)[2])
    names(label_sgl) <- colnames(X)
    label_sgl[colnames(cell_x)[which(sgl_results > 0)]] = "PACells"
    cat("Label Summary (B/P): ", table(label_sgl))
    sc_data <- seurat_object
    sc_data@meta.data$PACells_label <- factor(label_sgl, levels = c("PACells",
                                                                    "Background"))
    return(sc_data)
  }

  sc_data <-  identifyCellSub.Covar(X=simmtx, Y=phenotype,covariates =covariates, family = family,
                                    object = sc_dataset, candidatecell=cand_cell,
                                    cutoff = cutoff, resolution = res,
                                    dims_Neighbors = dims_Neighbors, dims_UMAP = dims_UMAP,
                                    group=group,
                                    batch = batch, Batch_variable = Batch_variable,
                                    dim_Batch = dim_Batch)

  return(sc_data)
}




#' PACells.RNA: identify clinical phenotype-associated cell states from scRNA-seq data using bulk phenotypes
#'
#'  \code{PACells.RNA} applies the PACells workflow to transcriptomic data. It links a phenotype defined on
#' bulk RNA-seq samples (e.g., case/control status, a continuous clinical variable, treatment response, or
#' survival outcome) to single-cell RNA-seq profiles and identifies a subset of cells most associated with
#' the bulk phenotype.
#'
#' The function expects two expression matrices (bulk and single-cell) and returns a Seurat object created
#' from the single-cell matrix with an additional metadata column PACells_label that labels each cell
#' as \code{"PACells"} or \code{"Background"}.
#'
#' \strong{Input requirements and alignment:}
#' \itemize{
#'   \item \code{sc_dataset} must be a gene-by-cell expression matrix (rows = genes, columns = cells).
#'   \item \code{bulk_dataset} must be a gene-by-sample expression matrix (rows = genes, columns = bulk samples).
#'   \item Row names (genes) are used to match bulk and single-cell data. Only shared genes are used.
#'   \item \code{phenotype} must be aligned to \code{colnames(bulk_dataset)} (same order).
#' }
#'
#' \strong{Workflow summary (high-level):}
#' \itemize{
#'   \item Construct a bulk-to-cell similarity matrix using the shared genes.
#'   \item Pre-screen a subset of candidate cells for efficiency and stability.
#'   \item Build a Seurat object from the single-cell matrix, perform standard RNA preprocessing and
#'   clustering/embedding (optionally with batch correction).
#'   \item Label phenotype-associated cells (\code{"PACells"}).
#' }
#'
#' @param sc_dataset A gene-by-cell expression matrix for single-cell RNA-seq data.
#' Rows are genes and columns are single cells. Row names should be gene symbols/IDs.
#' @param bulk_dataset A gene-by-sample expression matrix for bulk RNA-seq data.
#' Rows are genes and columns are bulk samples.
#' @param phenotype Phenotype annotation for bulk samples.
#' \itemize{
#'   \item For \code{family = "binomial"} and \code{family = "gaussian"}: a vector with length equal to the
#'   number of bulk samples, aligned to \code{colnames(bulk_dataset)}.
#'   \item For \code{family = "cox"}: a two-column matrix/data.frame with \code{time} and \code{status},
#'   aligned to bulk samples.
#' }
#' @param cutoff cutoff Maximum proportion of cells labeled as \code{"PACells"} among all single cells.
#' This parameter restricts the final number of selected cells. A smaller cutoff value (default \code{10\%}) is recommended
#' depending on the input data. Users can adjust it depending on dataset size and desired strictness.
#' @param screenRatio Pre-screen candidate cells that are strongly related to the phenotype.
#' This parameter is used to pre-screen the percentage of candidate cells based Bcor measure.
#' If a numeric value is provided, it directly specifies the proportion of cells to keep (default \code{0.80}).
#' An appropriate ratio (default \code{80\%}) is recommended depending on the input data.
#' If \code{screenRatio = "auto"}, PACells selects a suitable proportion by checking the stability of the
#' pre-screening results under repeated subsampling of bulk samples (controlled by \code{auto_subsample} and
#' \code{auto_stab_threshold}).
#' @param family Response type for the regression model. It depends on the type of the given phenotype and
#' can be \code{family = gaussian} for linear regression, \code{family = binomial} for classification,
#' or \code{family = cox} for Cox regression.
#' @param method Method for calculating the similarity matrix of bulk and single cells.
#' #' Options include \code{"KL"}, \code{"Pearson"}, and \code{"Spearman"}. The default is KL divergence.
#' @param level Logical value. If TRUE, calculating similarity based on an overall expression level profile.
#' If FALSE, calculating similarity based on relative expression patterns between features. The default is TRUE.
#' @param res the resolution parameter (default 1) in the FindClusters function, which is used to group cells.
#' @param dims_Neighbors Dimensions of reduction to construct nearest-neighbor graph.
#' @param dims_UMAP Dimensions of reduction to UMAP.
#' @param group Optional user-provided grouping vector for single cells. A vector of grouping information of
#' single cells, when is provided by the user, this grouping is used instead of clusters computed internally.
#'@param auto_subsample Subsampling proportion of bulk samples used when \code{screenRatio = "auto"}.
#' Default is \code{0.7}.
#'
#' @param auto_stab_threshold Stability threshold used when \code{screenRatio = "auto"}.
#' Default is \code{0.7}.
#'
#' @param batch Batch correction option applied to the scRNA-seq preprocessing step.
#' Supported options:
#' \itemize{
#'   \item \code{"none"}: no batch correction (default).
#'   \item \code{"harmony"}: Harmony integration applied after PCA.
#'   \item \code{"seurat"}: Seurat anchor-based integration (integrated assay followed by PCA).
#' }
#' When batch effects are evident (e.g., cells separate by sample/study in uncorrected embeddings or data are
#' combined across studies/platforms), we recommend using \code{batch = "harmony"}, which often provides a good
#' balance between mitigating batch-driven structure and preserving biological signal.
#'
#' @param Batch_variable Metadata column name in the internally created Seurat object indicating sample/batch IDs
#' used by \code{batch = "harmony"} or \code{batch = "seurat"}. Default is \code{"Sample"}.
#' If batch correction is used, users should ensure this variable is available for cells
#' externally and disable batch correction).
#'
#' @param dim_Batch Dimensions used for batch correction (Harmony/Seurat integration).
#'
#'
#' @return This function returns a Seurat-class object. It contains the results that cell state identified by PACells as strongly associated with the clinical phenotype.
#' A Seurat object created from \code{sc_dataset} with an additional metadata column \code{PACells_label}
#' (levels: \code{"PACells"} and \code{"Background"}). Cells labeled \code{"PACells"} are the subset identified as
#' most associated with the bulk phenotype.
#' @export

PACells.RNA <- function(sc_dataset, bulk_dataset, phenotype,
                    cutoff = 0.1, screenRatio = 0.8,
                    family = c("binomial", "gaussian", "cox")[1],
                    method = c("KL", "Pearson", "Spearman")[1],
                    level = TRUE,
                    res = 1, dims_Neighbors = 1:30, dims_UMAP = 1:30,
                    group = NULL, auto_subsample = 0.7,
                    auto_stab_threshold = 0.7,
                    batch = c("none", "harmony", "seurat")[1],
                    Batch_variable = "Sample",
                    dim_Batch = 1:30)
{

  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)

  if(method == "KL"){
    getSimilartyMatrix <- function(bulk_count, sc_count) {
      common <- intersect(rownames(bulk_count), rownames(sc_count))
      dataset0 <- cbind(bulk_count[common, ], sc_count[common,
      ])
      dataset1 <- preprocessCore::normalize.quantiles(dataset0,
                                                      keep.names = TRUE)
      bulk_norm <- dataset1[, 1:dim(bulk_count)[2]]
      sc_norm <- dataset1[, -c(1:dim(bulk_count)[2])]

      if (level){
        nbins <- seq(min(dataset1), max(dataset1), length.out = dim(dataset1)[1] *
                       2)
        sc_1 <- apply(sc_norm, 2, function(x) {
          hist(x, breaks = nbins, plot = FALSE)$counts
        })
        bulk_1 <- apply(bulk_norm, 2, function(x) {
          hist(x, breaks = nbins, plot = FALSE)$counts
        })
      } else {
        bulk_1 <- apply(bulk_norm, 2, robust_minmax)
        sc_1   <- apply(sc_norm,   2, robust_minmax)
      }
      bulk_1_list <- lapply(seq_len(ncol(bulk_1)), function(i) bulk_1[,
                                                                      i])
      KlDist <- function(x, y, base = 2) {
        x[x == 0] <- 1e-15
        y[y == 0] <- 1e-15
        x <- x/sum(x)
        y <- y/sum(y)
        D1 <- sum(x * log(x/y, base = base))
        D2 <- sum(y * log(y/x, base = base))
        D <- (D1 + D2)/2
        return(list(D1 = D1, D2 = D2, D = D))
      }
      KLMatrix1 <- t(as.data.frame(lapply(bulk_1_list, function(x) {
        apply(sc_1, 2, function(y) {
          KlDist(x, y)$D2
        })
      })))
      KL1Matrix <- max(KLMatrix1) - KLMatrix1
      rownames(KL1Matrix) <- colnames(bulk_count)
      colnames(KL1Matrix) <- colnames(sc_count)
      return(KL1Matrix)
    }
    simmtx <- getSimilartyMatrix(as.matrix(bulk_dataset), as.matrix(sc_dataset),
                                 level)

  }

  if(method == "Pearson"){

    getPearsonMatrix <- function(bulk_count, sc_count){

      dataset0 <- cbind(bulk_count,sc_count)
      dataset1 <- preprocessCore::normalize.quantiles(dataset0, keep.names = TRUE)
      bulk_count1 <- dataset1[,1:dim(bulk_count)[2]]
      sc_count1 <- dataset1[,-c(1:dim(bulk_count)[2])]
      pearsonMatrix <- cor(bulk_count1, sc_count1, method = c("pearson", "spearman")[1])

      return(pearsonMatrix)
    }

    simmtx <- getPearsonMatrix(bulk_count, sc_count)

  }

  if(method == "Spearman"){

    getSpearmanMatrix <- function(bulk_count, sc_count){

      dataset0 <- cbind(bulk_count,sc_count)
      dataset1 <- preprocessCore::normalize.quantiles(dataset0, keep.names = TRUE)

      bulk_count1 <- dataset1[,1:dim(bulk_count)[2]]
      sc_count1 <- dataset1[,-c(1:dim(bulk_count)[2])]

      spearsonMatrix <- cor(bulk_count1, sc_count1, method = c("pearson", "spearman")[2])

      return(spearsonMatrix)
    }

    simmtx <- getSpearmanMatrix(bulk_count, sc_count)

  }

  cat("Similarty Matrix Done","\n")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  getCandidateCell <- function(simmtx, phenotype, screenRatio = 0.8,
                               family = "binomial",
                               auto_ratios = c(0.7, 0.8, 0.9),
                               auto_B = 10, auto_subsample = 0.7,
                               auto_stab_threshold = 0.7){

    X <- simmtx
    Y <- phenotype

    jaccard <- function(a, b) {
      a <- unique(a); b <- unique(b)
      length(intersect(a, b)) / length(union(a, b))
    }

    run_one_screen <- function(Xb, Yb, ratio) {
      m <- dim(Xb)[2]
      filter_num <- max(1, (m * ratio))

      if(family == "binomial"){
        fit <- sbisis(X = Xb, Y = Yb,
                      candidate = filter_num,
                      method = "SBI-SIS-Pvalue",R=10)
      }
      if(family == "gaussian"){
        fit <- sbisis(X = Xb, Y = Yb,
                      candidate = filter_num,
                      method = "SBI-SIS-Pvalue",R=10)
      }
      if(family == "cox"){
        fit <- sbisis.surv(X = Xb, Y = Yb,
                           candidate = filter_num)
      }
      candidate_cell <- fit$candidate.set
    }

    autoPickRatio <- function(X, Y) {
      n <- nrow(X)

      stab_scores <- sapply(auto_ratios, function(r) {
        sets <- vector("list", auto_B)

        for (b in seq_len(auto_B)) {
          set.seed(10000 + b)
          idx <- sample(seq_len(n), size = max(2, floor(n * auto_subsample)))

          Xb <- X[idx, , drop = FALSE]
          if (family == "cox") {
            Yb <- Y[idx, , drop = FALSE]
          } else {
            Yb <- Y[idx]
          }

          sets[[b]] <- run_one_screen(Xb, Yb, r)
        }


        jac <- c()
        if (auto_B >= 2) {
          for (i in 1:(auto_B - 1)) {
            for (j in (i + 1):auto_B) {
              jac <- c(jac, jaccard(sets[[i]], sets[[j]]))
            }
          }
        }
        mean(jac, na.rm = TRUE)
      })

      ok <- which(stab_scores >= auto_stab_threshold)
      if (length(ok) > 0) {
        min(auto_ratios[ok])
      } else {
        0.8
      }
    }

    if (is.character(screenRatio) && screenRatio == "auto") {
      ratio <- autoPickRatio(X, Y)
    } else {
      ratio <- screenRatio
    }


    filter_num <- max(1, (dim(X)[2] * ratio))

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
  }
  cand_cell <- getCandidateCell(simmtx, phenotype, screenRatio=screenRatio,
                                family = family,
                                auto_subsample = auto_subsample,
                                auto_stab_threshold = auto_stab_threshold)

  cat("Pre-screening Done","\n")

  ##X is all all single cell similarty matrix,
  ##Y: phenotype label,
  ##counts: the matrix of all cell
  ##sc_peak: the peaks of single cell data
  ##candidatecell: the candidate cells by runing the SBISIS
  ##ratio: The highest ratio of cells were identified
  ##resolution: the parameter findclusters()
  identifyCellSub <- function (X, Y, family = c("binomial", "gaussian", "cox")[1], object,
                               candidatecell, cutoff = 0.1, resolution = 1,
                               dims_Neighbors = 1:30, dims_UMAP = 1:30, group = NULL,
                               batch = c("none", "harmony", "seurat")[1],
                               Batch_variable = "Sample",
                               dim_Batch = 1:30)
  {
    ratio <- cutoff
    seurat_object <- object
    seurat_object <- CreateSeuratObject(counts = as.matrix(seurat_object), assay = "RNA",
                                      min.cells = 1)
    seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
    seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = F)
    seurat_object <- ScaleData(seurat_object,verbose = F)
    seurat_object <- RunPCA(seurat_object, verbose = F)


    if(batch == "none") {
      reduction_use <- "pca"
    }
    if(batch == "harmony") {
      seurat_object <- RunHarmony(
        object = seurat_object, max_iter = 10,
        group.by.vars = Batch_variable,dims.use = dim_Batch,
        assay.use = 'RNA', project.dim = FALSE, verbose = F )
      reduction_use <- "harmony"
    }
    if (batch == "seurat") {
      obj_list <- SplitObject(seurat_object, split.by = Batch_variable)
      obj_list <- lapply(
        obj_list,
        function(x) {
          if ("RNA" %in% names(x@assays)) {
            DefaultAssay(x) <- "RNA"
          }
          x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
          x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = F)
          x <- ScaleData(x,verbose = F)
          x <- RunPCA(x, verbose = F)
          return(x)
        }
      )
      anchors <- FindIntegrationAnchors(
        object.list    = obj_list,
        reduction      = "rpca",
        dims           = dim_Batch
      )

      seurat_object <- IntegrateData(
        anchorset         = anchors,
        new.assay.name = "integrated",k.weight = 10,
        dims  =  dim_Batch
      )
      DefaultAssay(seurat_object) <- "integrated"
      seurat_object <- ScaleData(seurat_object, verbose = FALSE)
      seurat_object <- RunPCA(seurat_object, npcs = 100, verbose = FALSE)
      reduction_use <- "pca"

    }


    seurat_object <- FindNeighbors(seurat_object, reduction = reduction_use, dims = dims_Neighbors, verbose = F)
    seurat_object <- FindClusters(seurat_object,
                              resolution = resolution, verbose = F)
    seurat_object <- RunUMAP(seurat_object, reduction = reduction_use, dims = dims_UMAP, verbose = F)

    if(!is.null(group)){

      index <- group[candidatecell]

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

  sc_data <-  identifyCellSub(X=simmtx, Y=phenotype, family = family,
                              object = sc_dataset, candidatecell=cand_cell,
                              cutoff = cutoff, resolution = res,
                              dims_Neighbors = dims_Neighbors, dims_UMAP = dims_UMAP,
                              group = NULL, batch = batch, Batch_variable = Batch_variable,
                              dim_Batch = dim_Batch
                              )

  return(sc_data)
}

#' Robust min-max normalization for a numeric vector. This helper function rescales a numeric vector to
#' approximately the [0, 1] range using min-max normalization.
#' @param x A numeric vector to be rescaled.
#' @param eps A small positive constant added to the output to avoid exact zeros
#' (default 1e-8).
#' @return A numeric vector of the same length as x, rescaled by min-max normalization
#' and shifted by eps
#' @noRd
robust_minmax <- function(x, eps = 1e-8) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (!is.finite(rng) || rng < 1e-15) {
    y <- rep(1, length(x))
  } else {
    y <- (x - min(x, na.rm = TRUE)) / rng
  }
  y + eps
}


#These screening codes is obtained from https://doi.org/10.1080/01621459.2018.1462709.
##  Pre-screening candidate cells based on Bcor.
## In PACells, this function is used to pre-screen a subset of candidate single cells
## from a bulk–cell similarity matrix.
##
##' This function ranks features (e.g., candidate cells) in \code{X} by their marginal association
##' with the phenotype matrix/vector \code{Y} using Bcor (SBI/BI). It then returns
##' the top \code{candidate} features as a candidate set.
##'
##' Typical usage in PACells:
##' \itemize{
##'   \item \code{X}: a similarity matrix with rows = bulk samples and columns = single cells.
##'   \item \code{Y}: bulk phenotype aligned to the same bulk samples (rows of \code{X}).
##' }
##' @param Y Y A numeric vector or matrix of phenotypes for bulk samples. If a matrix, rows correspond
##'   to samples and columns correspond to phenotype variables.
##' @param X X A numeric matrix (e.g., similarity matrix). Rows correspond to bulk samples and columns
##'   correspond to features to be screened (e.g., single cells).
##' @param candidate candidate Candidate set size. Can be:
##'   \itemize{
##'     \item a numeric value specifying the number of candidates to keep;
##'     \item \code{"small"}: keeps approximately \code{round(n/log(n))} candidates (n = number of samples);
##'     \item \code{"large"}: keeps a large candidate set (up to all available features).
##'   }
##' @param method Method for screening procedure.
##' @param parms parameters list only available when method = "lm" or "gam".
##' It contains three parameters: d1, d2, and df.
##' d1 is the number of initially selected variables,
##' d2 is the number of variables collection size added in each iteration.
##' df is degree freedom of basis in generalized additive models playing a role only when method = "gam".
##' Default:  parms = list(d1 = 5, d2 = 5, df = 3)
##' @param R Number of permutations used.
##' @param parallel Logical; whether to enable parallel screening using \code{snowfall}.
##' @param ncore Integer; number of CPU cores used when \code{parallel = TRUE}. If NULL, snowfall uses its default.
##' @return A list with:
##' \itemize{
##'   \item \code{candidate.set}: indices of selected candidate features (columns of \code{X});
##'   \item \code{candidate.data}: a list containing \code{X} restricted to selected features and the (possibly updated) \code{Y};
##'   \item \code{candidate.size}: number of selected candidates.
##' }
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
      Xlastpickout <- ids[chooseids]

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


#' BI: Bcor-based statistic
#'
#' This function is obtained from https://doi.org/10.1080/01621459.2018.1462709
#' and included here as part of the PACells
#' implementation.
#'
#' @param x A numeric vector or matrix. If a matrix, rows correspond to samples and columns correspond
#' to cells. In PACells, \code{x} typically represents a candidate cell across bulk samples.
#' @param y A numeric vector or matrix. If a matrix, rows correspond to samples and columns correspond
#' to phenotype variables. In PACells, \code{y} represents the bulk phenotype aligned to samples.
#' @param R Number of permutations. Use \code{R = 0} to compute the statistic only; use \code{R > 0}
#' to enable permutation-based evaluation.
#' @param seed Integer random seed used.
#' @noRd
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


#' This function is obtained from the Bcor/SBI-SIS framework and included in PACells .
#' @param x a numeric vector.
#' @param y a numeric vector  (aligned to \code{x} by sample order).
#' @param R Number of permutations (0 means no permutation).
#' @param seed Random seed used.
#' @return A numeric value returned by the compiled UBI routine.
#' @noRd
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
#' This function is adapted from the Bcor/SBI-SIS framework
#' and used in PACells as a basic association score between x and y. It relies on BI() for
#' the general case, and UBI() for the fast univariate case.
#' @aliases SBI
#' @param x a numeric matirx
#' @param y a numeric matirx (aligned to \code{x} by sample order).
#' @param fast Logical; if TRUE, use the univariate implementation (\code{UBI}) for speed.
#' @return A numeric SBI score.
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

##' sbisis.surv: pre-screening for survival outcomes
##'
##' ##' This function is adapted from the Bcor/SBI-SIS framework
##' and used in PACells to pre-screen candidate features (e.g., cells) when the bulk phenotype is survival data.
##' It ranks columns of \code{X} by a survival-version screening score and returns the top candidates.
##' @param Y a numeric matirx(first column should be event time, second column should be survival status) or Surv object
##' @param X A numeric matrix with rows = samples and columns = features to be screened.
##' @param candidate size of candidate set.
##' @param standized standized Logical; whether to standardize covariates in \code{X} before screening.
##' @return A list containing \code{candidate.set}, \code{candidate.data}, and \code{candidate.size}.
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

# This function is adapted from the Bcor/SBI-SIS framework
# and included in PACells as a low-level routine for survival screening.
# @param x Ordered covariate vector
# @param t ordered survival event time
# @param delta ordered survival event status
# @param Sc Survival function values from a Kaplan–Meier fit (expanded to match t if needed).
# @param n Sample size
# @useDynLib PACells
# @export
# @noRd
#
SRCT <- function(x, t, delta, Sc, n){
  RCT <- numeric(1)
  RCT<-.C("SRCT", as.double(t(x)), as.double(t(t)), as.double(t(delta)),
          as.double(t(Sc)), as.integer(n), RC=as.double(RCT))
  return(RCT$RC)
}

#' bcor.surv: survival ball-correlation score
#'
#' This function is adapted from the Bcor/SBI-SIS framework
#' and provides a convenient wrapper to compute the survival association score between a covariate \code{x}
#' and survival outcome (\code{time}, \code{status}). It standardizes \code{x}, constructs the survival function,
#' and calls \code{SRCT()}.
#' @param time Survival/event time.
#' @param status Event status.
#' @param x Numeric vector (aligned to \code{time} and \code{status}).
#' @return Survival Ball Correlation. A numeric survival association score.
#' @import stats
#' @import utils
#' @import survival
#' @export
#' @noRd
#' @examples
#' data(survdat)
#' bcor.surv(time = survdat[["time"]], status = survdat[["status"]], x = survdat[,3])
bcor.surv <- function(time, status, x) {
  n = length(time)
  ord.t = sort(time)
  ix = order(time)
  ord.delta = status[ix]
  x = scale(x[ix])

  fitc = survfit(Surv(time,1-status)~1)
  Sc = fitc$surv
  if(length(unique(ord.t)) != n) {
    rep_num = as.data.frame(table(ord.t))[, "Freq"]
    Sc = mapply(function(x, y) {
      rep(x, y)
    }, Sc, rep_num, SIMPLIFY = FALSE)
    Sc= unlist(Sc)
  }

  SRCT(x = x,t = ord.t, delta = ord.delta, Sc = Sc, n = n)
}

