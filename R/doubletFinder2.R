
doubletFinder2 <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
          sct = FALSE, annotations = NULL, returndoublets=TRUE, seu_wdoublets=NULL, 
          maxN=50, seed=1234)  {
  require(Seurat)
  require(fields)
  require(KernSmooth)
  set.seed(seed)
  is.new.run <- is.null(seu_wdoublets)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts[, real.cells]
    n_real.cells <- length(real.cells)
    
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      
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
      colnames(doublets) <- paste0(as.character(doublet_types1), 
                                   "__", 
                                   as.character(doublet_types2),
                                   "XX", 1:nrow(db))
      
    } else {
      n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
      print(paste("Creating", n_doublets, "artificial doublets...", 
                  sep = " "))
      real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
      real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
      doublets <- (data[, real.cells1] + data[, real.cells2])/2
      colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    }
    data_wdoublets <- cbind(data, doublets)
    
    orig.commands <- seu@commands
    if (sct == FALSE & is.null(seu_wdoublets)) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
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
      seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                              npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                              verbose = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      if(!returndoublets) rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE & is.null(seu_wdoublets)) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      if(!returndoublets) rm(seu_wdoublets)
      gc()
    }
    if(is.new.run & !is.null(annotations)){
      seu_wdoublets$dubid <-  'NA'
      doubletidx <- grep("__.*XX", Cells(seu_wdoublets))
      seu_wdoublets$dubid[doubletidx] <- gsub("XX[0-9]*$", "", Cells(seu_wdoublets)[doubletidx])
    }
    if(!is.null(seu_wdoublets)){
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      if(!returndoublets) rm(seu_wdoublets)
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol =2))
      colnames(neighbor_types) <- c('homotypic', 'heterotypic')
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if (!is.null(annotations)) {
        for (ct in unique(as.character(annotations))) {
          neighbors_that_are_doublets = neighbors[neighbors >  n_real.cells]
          if (length(neighbors_that_are_doublets) > 0) {
            idx <- neighbors_that_are_doublets -  n_real.cells
            neighbor_types[i, ] <- table(
              factor(as.character(doublet_types1[idx]) == as.character(doublet_types2[idx]),
                     levels=c('TRUE', 'FALSE')))
            neighbor_types[i, ] <- neighbor_types[i, ]/sum(neighbor_types[i, ])
          }
          else {
            neighbor_types[i, ] <- NA
          }
        }
      }
    }  
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data),  1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    if (!is.null(annotations)) {
      # colnames(neighbor_types) = levels(doublet_types1)
      for (ct in colnames(neighbor_types)) {
        seu@meta.data[, paste("DF.doublet.contributors", 
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[,ct]
      }
    }
    if(returndoublets){
      return(list("seu"=seu, 
                  "seu_wdoublets"=seu_wdoublets,
                  "pca.coord"=pca.coord))
    } else {
      return(seu)
    }
  }
}
