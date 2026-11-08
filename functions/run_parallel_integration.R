#!/usr/bin/env Rscript
# Independent stages share the same HVGs, PCA, cell order and dataset batch model.
suppressPackageStartupMessages({library(Seurat); library(jsonlite)})
a <- commandArgs(TRUE); stopifnot(length(a) == 2L)
run <- normalizePath(a[1]); stage <- a[2]
stopifnot(stage %in% c('pca','rpca','harmony','raw_umap','collect','cluster'))
Sys.setenv(BRAINOMICS_RUN_ROOT=run)
script_file <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1L])
repo <- Sys.getenv('BRAINOMICS_REPO_ROOT', unset=dirname(dirname(normalizePath(script_file))))
setwd(repo)
source('functions/data_paths.R'); source('functions/dataset_metadata.R')
source('functions/processed_object.R'); source('functions/integration.R')
configure_complete_integration_future()
seed <- as.integer(Sys.getenv('BRAINOMICS_INTEGRATION_SEED','20260730'))
dims <- seq_len(as.integer(Sys.getenv('BRAINOMICS_INTEGRATION_DIMS','50')))
n_features <- as.integer(Sys.getenv('BRAINOMICS_VARIABLE_FEATURES','3000'))
out <- brainomics_results_dir(); cp <- file.path(out,'r_checkpoints')
dir.create(cp,recursive=TRUE,showWarnings=FALSE)
policy_hash <- processed_file_sha256(file.path(run,'named_gene_panel_policy.json'))
input <- file.path(cp,'pca_input.rds'); pca_file <- file.path(cp,'objects_pca.rds')
result_file <- file.path(cp,paste0(stage,'_reductions.rds'))
contract <- list(Seed=seed, Dimensions=dims, Variable_Features=n_features,
  Batch_Model='dataset', Reference_Datasets=reference_datasets(),
  Panel_Policy_SHA256=policy_hash)
atomic_rds <- function(value, path) {
  stopifnot(!file.exists(path))
  temp <- paste0(path,'.tmp.',Sys.getpid())
  saveRDS(value,temp,compress=FALSE); stopifnot(file.rename(temp,path))
}
timing <- function(event) {
  write.table(data.frame(Stage=stage,Event=event,
    Timestamp_UTC=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC')),
    file.path(out,paste0('timing_',stage,'.tsv')),sep='\t',quote=FALSE,
    row.names=FALSE,col.names=identical(event,'start'),append=!identical(event,'start'))
}
check_common <- function(common) {
  stopifnot(identical(common$contract,contract), nrow(common$pca)==2602031,
    ncol(common$pca)==length(dims), !anyDuplicated(rownames(common$pca)),
    identical(rownames(common$pca),rownames(common$metadata)),
    length(common$features)==n_features,
    identical(sort(unique(as.character(common$metadata$Dataset))),sort(reference_datasets())),
    identical(as.character(common$metadata$Integration_Batch_ID),
              paste0('dataset:',common$metadata$Dataset)))
}
make_umap <- function(embedding, key) {
  set.seed(seed)
  RunUMAP(embedding[,dims,drop=FALSE], assay='RNA', reduction.key=key,
          seed.use=seed, verbose=FALSE)
}
timing('start'); set.seed(seed)
if(stage == 'cluster') {
  future::plan('sequential')
  common <- readRDS(input); check_common(common)
  result <- readRDS(file.path(cp,'rpca_reductions.rds'))
  stopifnot(identical(result$contract,contract))
  embedding <- Embeddings(result$reductions$integrated.rpca)
  cells <- rownames(embedding)
  stopifnot(identical(dim(embedding),c(2602031L,50L)),
    identical(cells,rownames(common$pca)), !anyDuplicated(cells), all(is.finite(embedding)))
  outputs <- file.path(out,c('rpca_knn_neighbor.rds','rpca_symmetric_weighted_knn_graph.rds',
    'cluster_assignments.rds','cluster_sizes.tsv','clustering_contract.json'))
  stopifnot(!any(file.exists(outputs)))
  rm(common,result); gc()
  k <- 20L; resolution <- 2
  set.seed(seed)
  neighbor <- Seurat::FindNeighbors(embedding,k.param=k,return.neighbor=TRUE,
    compute.SNN=FALSE,nn.method='annoy',annoy.metric='euclidean',n.trees=50L,verbose=TRUE)
  stopifnot(identical(dim(neighbor@nn.idx),c(length(cells),k)),
    identical(neighbor@cell.names,cells),all(is.finite(neighbor@nn.dist)))
  atomic_rds(neighbor,outputs[1L])
  rm(embedding); gc()
  idx <- neighbor@nn.idx; distance <- neighbor@nn.dist
  scale <- pmax(distance[,k],.Machine$double.eps)
  src <- rep(seq_len(nrow(idx)),each=k); dst <- as.integer(t(idx))
  d <- as.numeric(t(distance)); s <- rep(scale,each=k); keep <- src!=dst
  graph <- Matrix::sparseMatrix(i=src[keep],j=dst[keep],
    x=exp(-((d[keep]/s[keep])^2)),dims=c(length(cells),length(cells)),dimnames=list(cells,cells))
  rm(neighbor,idx,distance,scale,src,dst,d,s,keep); gc()
  graph <- as(Matrix::drop0(graph+Matrix::t(graph)),'Graph')
  stopifnot(isSymmetric(graph),all(is.finite(graph@x)),all(graph@x>0))
  atomic_rds(graph,outputs[2L])
  set.seed(seed)
  partition <- FindClusters(graph,algorithm=1L,resolution=resolution,n.start=10L,
    n.iter=10L,random.seed=seed,group.singletons=TRUE,verbose=TRUE)
  labels <- paste0('C',as.character(partition[[1L]]))
  stopifnot(length(labels)==length(cells),!anyNA(labels))
  atomic_rds(data.frame(Cell=cells,Cluster=labels),outputs[3L])
  write.table(as.data.frame(table(Cluster=labels)),outputs[4L],sep='\t',quote=FALSE,row.names=FALSE)
  write_json(list(state='PASS',cells=length(cells),clusters=length(unique(labels)),
    dimensions=50,k_neighbors=k,n_trees=50,metric='euclidean',
    graph='symmetric adaptive Gaussian weighted kNN; reciprocal weights summed',
    algorithm='Louvain (Seurat algorithm 1)',resolution=resolution,n_start=10,n_iter=10,
    seed=seed,old_assignments_reused=FALSE,biological_annotation=FALSE,contract=contract),
    outputs[5L],pretty=TRUE,auto_unbox=TRUE)
} else if(stage == 'pca') {
  stopifnot(!file.exists(input), !file.exists(pca_file))
  object <- load_processed_object(file.path(cp,'objects_hvg.rds'))
  hvg <- object@misc$BrainOmics_Integration_Checkpoint
  stopifnot(hvg$Stage=='hvg', hvg$Seed==seed, hvg$Variable_Features==n_features,
    identical(as.integer(hvg$Dimensions),dims), hvg$Batch_Model=='dataset',
    identical(hvg$Reference_Datasets,reference_datasets()),
    ncol(object)==2602031, nrow(object)==24659,
    length(VariableFeatures(object))==n_features)
  object <- ScaleData(object,features=VariableFeatures(object),verbose=FALSE)
  object <- RunPCA(object,features=VariableFeatures(object),npcs=max(dims),
                    seed.use=seed,verbose=FALSE)
  common <- list(pca=Embeddings(object,'pca'),
    metadata=object[[]][,c('Dataset','Integration_Batch_ID'),drop=FALSE],
    features=VariableFeatures(object), contract=contract)
  check_common(common)
  object@misc$Parallel_Integration_Contract <- contract
  save_processed_object(object,pca_file,compression='gzip',validate_reload=FALSE)
  atomic_rds(common,input)
} else {
  common <- readRDS(input); check_common(common)
  input_sha <- processed_file_sha256(input)
  if(stage == 'collect') {
    object <- load_processed_object(pca_file)
    stopifnot(identical(object@misc$Parallel_Integration_Contract,contract),
      identical(Embeddings(object,'pca'),common$pca))
    versions <- list()
    for(method in c('rpca','harmony','raw_umap')) {
      result <- readRDS(file.path(cp,paste0(method,'_reductions.rds')))
      stopifnot(identical(result$contract,contract),result$input_sha256==input_sha)
      for(name in names(result$reductions)) {
        reduction <- result$reductions[[name]]
        stopifnot(identical(rownames(Embeddings(reduction)),colnames(object)),
                  all(is.finite(Embeddings(reduction))))
        object[[name]] <- reduction
      }
      versions[[method]] <- result$versions
    }
    identity <- integration_method_identity()
    identity <- identity[identity$Display_Name != 'scVI',,drop=FALSE]
    stopifnot(all(identity$Reduction %in% names(object@reductions)))
    write_tsv(identity,file.path(out,'integration_method_identity_r.tsv'))
    write_json(versions,file.path(out,'parallel_integration_versions.json'),
               pretty=TRUE,auto_unbox=TRUE)
    # All R integration is complete; final objects retain counts/data and reductions.
    # The reproducible 62 GB global scaled layer stays in the PCA checkpoint only.
    object[['RNA']]$scale.data <- NULL
    save_processed_object(object,file.path(out,'objects_integrated_r.rds'),
                           compression='gzip',validate_reload=FALSE)
  } else {
    stopifnot(!file.exists(result_file))
    packages <- c('Seurat','SeuratObject',if(stage=='harmony') 'harmony')
    versions <- setNames(lapply(packages,function(p) as.character(packageVersion(p))),packages)
    latent_file <- file.path(cp,paste0(stage,'_latent.rds'))
    if(stage != 'raw_umap' && file.exists(latent_file)) {
      cached <- readRDS(latent_file)
      stopifnot(identical(cached$contract,contract),cached$input_sha256==input_sha,
                identical(cached$versions,versions))
      reduction <- cached$reduction
      timing('latent_resume')
    } else if(stage == 'rpca') {
      object <- load_processed_object(pca_file)
      stopifnot(identical(object@misc$Parallel_Integration_Contract,contract),
                identical(Embeddings(object,'pca'),common$pca))
      timing('latent_start')
      object <- IntegrateLayers(object=object,method=RPCAIntegration,
        orig.reduction='pca',new.reduction='integrated.rpca',verbose=FALSE)
      reduction <- object[['integrated.rpca']]
      rm(object); gc(FALSE)
      timing('latent_end')
      atomic_rds(list(reduction=reduction,contract=contract,input_sha256=input_sha,
                      versions=versions),latent_file)
    } else if(stage == 'harmony') {
      stopifnot(as.character(utils::packageVersion('harmony')) == '2.0.5')
      cores <- as.integer(Sys.getenv('BRAINOMICS_HARMONY_CORES','8'))
      timing('latent_start')
      embedding <- harmony::RunHarmony(data_mat=common$pca,
        meta_data=common$metadata,vars_use='Integration_Batch_ID',
        theta=NULL,lambda=NULL,sigma=0.1,nclust=NULL,max_iter=10,
        ncores=cores,return_object=FALSE,verbose=TRUE,
        .options=harmony::harmony_options(tau=0,block.size=0.05,
          max.iter.cluster=20,epsilon.cluster=1e-5,epsilon.harmony=0.01))
      rownames(embedding) <- rownames(common$pca)
      colnames(embedding) <- paste0('harmony_',seq_len(ncol(embedding)))
      reduction <- CreateDimReducObject(embeddings=embedding,key='harmony_',assay='RNA')
      timing('latent_end')
      atomic_rds(list(reduction=reduction,contract=contract,input_sha256=input_sha,
                      versions=versions),latent_file)
    }
    timing('umap_start')
    reductions <- switch(stage,
      rpca=list(integrated.rpca=reduction,umap.rpca=make_umap(Embeddings(reduction),'UMAP_')),
      harmony=list(integrated.harmony=reduction,umap.harmony=make_umap(Embeddings(reduction),'UMAP_')),
      raw_umap=list(umap.unintegrated=make_umap(common$pca,'UMAP_')))
    timing('umap_end')
    stopifnot(all(vapply(reductions,function(x)
      identical(rownames(Embeddings(x)),rownames(common$pca)) &&
      all(is.finite(Embeddings(x))),logical(1))))
    atomic_rds(list(reductions=reductions,contract=contract,input_sha256=input_sha,
                    versions=versions),result_file)
  }
}
timing('end'); message(Sys.time(),' completed ',stage)
