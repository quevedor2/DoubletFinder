
.createDoublets <- function(data, annotations=NULL, maxN=NULL, pN=0.25){
  if(is.null(annotations)){
    real.cells <- colnames(data)
    n_real.cells <- length(real.cells)
    
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste(sample(real.cells, size=length(real.cells1), replace=T),
                                "..__doubletXX", 1:n_doublets, sep = "")
    
  } else {
    anno_spl <- lapply(split(names(annotations), annotations), function(i){
      sample(i, (pN * length(i)))
    })
    n <- if(is.null(maxN)) quantile(sapply(anno_spl, length), 0.33) else maxN
    uniq_combs <- combn(seq_along(anno_spl), m=2)
    db <- apply(uniq_combs, 2, function(ij){
      i1 <- anno_spl[[ij[1]]]
      i2 <- anno_spl[[ij[2]]]
      data.frame("rc1"=sample(i1, n, replace=T),
                 "rc2"=sample(i2, n, replace=T)) %>%
        dplyr::mutate("dt1"=as.character(annotations[rc1]),
                      "dt2"=as.character(annotations[rc2])) %>%
        dplyr::filter(dt1!=dt2)
    }) %>% do.call(rbind, .)
    print(paste("Creating", nrow(db), "artificial HETEROTYPIC doublets...", 
                sep = " "))
    real.cells1 <- db$rc1
    real.cells2 <- db$rc2
    doublet_types1 <- db$dt1
    doublet_types2 <- db$dt2
    
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste0(sample(names(annotations), size=length(doublet_types1), replace=T),
                                 "..",
                                 as.character(doublet_types1), 
                                 "__", 
                                 as.character(doublet_types2),
                                 "XX", 1:nrow(db))
  }
  return(doublets)
}

.preprocessSeu <- function(data_wdoublets, sct, PCs, resolution, annotations=NULL, orig.commands=NULL){
  print("Creating Seurat object...")
  seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
  
  if (sct == TRUE) {
    require(sctransform)
    print("Running SCTransform...")
    seu_wdoublets <- SCTransform(seu_wdoublets)
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
  } else if (sct == FALSE){
    print("Normalizing Seurat object...")
    seu_wdoublets <- NormalizeData(seu_wdoublets, 
                                   normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                   scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                   margin = orig.commands$NormalizeData.RNA@params$margin)
    print("Finding variable genes...")
    seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                          selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                          loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                          clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                          mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                          dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                          num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                          binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                          nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                          mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                          dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
    print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                               model.use = orig.commands$ScaleData.RNA$model.use, 
                               do.scale = orig.commands$ScaleData.RNA$do.scale, 
                               do.center = orig.commands$ScaleData.RNA$do.center, 
                               scale.max = orig.commands$ScaleData.RNA$scale.max, 
                               block.size = orig.commands$ScaleData.RNA$block.size, 
                               min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, 
                            features = orig.commands$ScaleData.RNA$features, 
                            npcs = length(PCs), 
                            rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                            verbose = FALSE)
  }
  
  if(!is.null(annotations)){
    seu_wdoublets$annotations <- 'doublets'
    seu_wdoublets@meta.data[names(annotations), 'annotations'] <- annotations
  }
  
  
  print(paste0("Finding small resolution clusters: resolution = ", resolution))
  seu_wdoublets <- FindNeighbors(seu_wdoublets, reduction = "pca", dims = 1:length(PCs)) %>%
    FindClusters(., resolution = resolution,  cluster.name = "small_clus") %>%
    RunUMAP(., reduction = "pca", dims = 1:length(PCs), n.neighbors=30L,
            min.dist = 0.1)
  
  print("Marking generated doublets...")
  seu_wdoublets$dubid <-  'NA'
  doubletidx <- grep("__.*XX", Cells(seu_wdoublets))
  seu_wdoublets$dubid[doubletidx] <- gsub("XX[0-9]*$", "", Cells(seu_wdoublets)[doubletidx]) %>%
    gsub("^.*\\.\\.", "", .)
  return(seu_wdoublets)
}


doubletFinder2 <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
          sct = FALSE, annotations = NULL, seu_wdoublets=NULL, 
          maxN=50, seed=1234)  {
  require(Seurat)
  require(fields)
  require(KernSmooth)
  set.seed(seed)
  is.new.run <- is.null(seu_wdoublets)
  if (reuse.pANN) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  

  ## Extract data and real-cell barcodes
  real.cells <- rownames(seu@meta.data)
  data <- seu@assays$RNA$counts[, real.cells]
  n_real.cells <- length(real.cells)
  
  ## Validation and formatting annotation vector
  if (!is.null(annotations)) {
    stopifnot(typeof(annotations) == "character")
    stopifnot(length(annotations) == length(Cells(seu)))
    stopifnot(!any(is.na(annotations)))
    annotations <- factor(annotations)
    names(annotations) <- Cells(seu)
  }
  
  ## Create doublets with or without annotations
  doublets <- .createDoublets(data, annotations =annotations, maxN=maxN, pN=pN)
  data_wdoublets <- cbind(data, doublets)
  orig.commands <- seu@commands
  
  ## Preprocess seurat object if not given
  if(is.null(seu_wdoublets)){
    seu_wdoublets <- .preprocessSeu(data_wdoublets, sct=sct, PCs=PCs, 
                                    resolution=resolution, annotations=annotations,
                                    orig.commands=orig.commands)
  }
  pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,PCs]  ## PCA matrix
  cell.names <- rownames(seu_wdoublets@meta.data)
  
  ## Calculate a distance matrix based on PC-coordinates
  print("Calculating PC distance matrix...")
  dist.mat <- fields::rdist(pca.coord) #Distance matrix of cells + doublets
  colnames(dist.mat) <- rownames(dist.mat) <- rownames(pca.coord)
  
  ## Compute counts of nearest neighbours 
  print("Computing pANN...")
  pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                               ncol = length(pK)))
  rownames(pANN) <- real.cells
  colnames(pANN) <- paste0("pANN_", pK)
  k <- round(n_real.cells * pK)

  for (i in 1:n_real.cells) {
    neighbors <- order(dist.mat[, i])
    pann_k <- sapply(k, function(ki){
      neighbors2 <- neighbors[2:(ki + 1)]
      length(which(neighbors2 > n_real.cells))/ki
    })
    pANN[i,] <- pann_k
  }  
  
  print("Classifying doublets..")
  classifications <- as.data.frame(matrix("Singlet", nrow = n_real.cells, 
                               ncol = length(pK)))
  rownames(classifications) <- real.cells
  colnames(classifications) <- paste0("pANN_", pK)
  for (i in seq_along(k)) {
    classifications[order(pANN[1:n_real.cells,i], decreasing = TRUE)[1:nExp], i] <- "Doublet"
  }
  
  
  seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data),  1]
  seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
  
  
  return(list("seu_wdoublets"=seu_wdoublets,
              "pca.coord"=pca.coord,
              'pANN'=pANN,
              'classifications'=classifications))
}
