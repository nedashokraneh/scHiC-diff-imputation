source("R/coolR.R")
source("R/utilities.R")
library('rhdf5')
library(pbapply)
library(data.table)
library(InteractionSet)

options(scipen = 999)

GI2tab <- function(GI, chrom){
  anchor1_GRanges = ranges(anchors(GI)$first)
  anchor2_GRanges = ranges(anchors(GI)$second)
  anchor1_starts = start(anchor1_GRanges) - 1
  anchor2_starts = start(anchor2_GRanges) - 1
  tab = data.frame(chr = chrom,
                   region1 = anchor1_starts,
                   region2 = anchor2_starts,
                   IF = mcols(GI)$count)
  tab
}

GIs2ISet <- function(GIs, chrom, resolution, max_distance){
  tabs = lapply(seq(length(GIs)), function(i) GI2tab(GIs[[i]], chrom))
  print(colnames(tabs[[1]]))
  for (i in c(1:length(tabs))){
    names(tabs[[i]])[names(tabs[[i]]) == 'IF'] = paste0('IF', i)
    tabs[[i]] = tabs[[i]][tabs[[i]]$region2 != tabs[[i]]$region1,]
    if (!is.null(max_distance)){
      tabs[[i]] = tabs[[i]][tabs[[i]]$region2 - tabs[[i]]$region1 <= max_distance,]
    }
  }
  dts = lapply(tabs, function(tab) data.table(tab, key = c('chr', 'region1', 'region2')))
  merged_dt = Reduce(function(dt1,dt2) {
    merge.data.table(dt1, dt2, all = TRUE)},
    dts)
  #merged_dt = merged_dt[merged_dt$region1 != merged_dt$region2]
  #if (!is.null(max_distance)){
  #  merged_dt = merged_dt[(merged_dt$region2 - merged_dt$region1) <= max_distance]
  #}
  source_gr = GRanges(seqnames = merged_dt$chr, 
                      ranges = IRanges(merged_dt$region1, 
                                       end = merged_dt$region1+resolution))
  target_gr = GRanges(seqnames = merged_dt$chr,
                      ranges = IRanges(merged_dt$region2, 
                                       end = merged_dt$region2+resolution))
  GI = GInteractions(source_gr, target_gr)
  counts = as(merged_dt[,4:ncol(merged_dt)], 'matrix')
  counts[is.na(counts)] = 0
  #colnames(counts) <- seq_len(ncol(counts))
  iset = InteractionSet(list(counts = counts), GI)
  return(iset)
}

read_tab_hic <- function(tab_path, file_id, colnames, chromosomes){
  dt = fread(tab_path, col.names = colnames)
  dt = dt[dt$chr1 == dt$chr2]
  dt$chr2 = NULL
  setnames(dt, 'chr1', 'chr')
  setnames(dt, 'count', paste0('count', file_id))
  dt = dt[dt$chr %in% chromosomes]
  setkeyv(dt, c('chr', 'pos1', 'pos2'))
  dt
}

load_sparse_iset_from_mcool <- function(mcool_paths, chroms, 
                                        resolution, chrom_name, 
                                        max_distance, mcool){
  if (mcool){
    GIs = pblapply(c(1:length(mcool_paths)),
                   function(i) read.cool(mcool_paths[[i]], res = resolution, 
                                         chr1 = chroms[[i]], chr2 = chroms[[i]]))
  }
  else{
    GIs = pblapply(c(1:length(mcool_paths)),
                   function(i) read.cool(mcool_paths[[i]], res = NULL, 
                                         chr1 = chroms[[i]], chr2 = chroms[[i]]))
  }
  iset = GIs2ISet(GIs, chrom_name, resolution, max_distance)
  iset
}

load_sparse_iset_from_tab <- function(tab_paths, colnames, chromosomes, resolution){
  dts = pblapply(c(1:length(tab_paths)), function(file_id) read_tab_hic(tab_paths[[file_id]], file_id,
                                                                        colnames, chromosomes))
  merged_dt =  Reduce(function(dt1,dt2){merge.data.table(dt1,dt2,all = TRUE)}, 
                      dts)
  merged_dt[is.na(merged_dt)] = 0
  source_gr = GRanges(seqnames = merged_dt$chr, 
                      ranges = IRanges(merged_dt$pos1, end = merged_dt$pos1 + resolution))
  target_gr = GRanges(seqnames = merged_dt$chr,
                      ranges = IRanges(merged_dt$pos2, end = merged_dt$pos2 + resolution))
  GI = GInteractions(source_gr, target_gr)
  iset = InteractionSet(as(merged_dt[,4:ncol(merged_dt)], 'matrix'), GI)
  elementMetadata(iset) = merged_dt$chr
  iset
}

sparse2dense_iset <- function(sparse_iset){
  rg = regions(sparse_iset)
  num_samples = dim(sparse_iset)[2]
  cm_list = lapply(c(1:num_samples), function(i) inflate(sparse_iset, rg, rg, assay=1, sample = i))
  cm_list = lapply(cm_list, function(cm){
    cm_matrix = as.matrix(cm)
    cm_matrix[is.na(cm_matrix)] = 0
    as.matrix(cm) = cm_matrix
    cm
  })
  dense_isets = lapply(cm_list, function(cm) deflate(cm, use.zero = TRUE))
  dense_iset = do.call(cbind, dense_isets)
  colnames(dense_iset) = colnames(sparse_iset)
  dense_iset
}

# faster when number of cells grows for single-cell data

load_dense_iset_from_mcool <- function(mcool_paths, chroms1, chroms2, 
                                       resolution, mcool){
  if (mcool){
    GIs = pblapply(c(1:length(mcool_paths)),
                   function(i) read.cool(mcool_paths[[i]], res = resolution, 
                                         chr1 = chroms1[[i]], chr2 = chroms2[[i]]))
  }
  else{
    GIs = pblapply(c(1:length(mcool_paths)),
                   function(i) read.cool(mcool_paths[[i]], res = NULL, 
                                         chr1 = chroms1[[i]], chr2 = chroms2[[i]]))
  }
  rg = regions(GIs[[1]])
  dense_isets = pblapply(GIs, function(GI) {
    cm = inflate(GI, rg, rg, fill = mcols(GI)$count)
    cm_matrix = as.matrix(cm)
    cm_matrix[is.na(cm_matrix)] = 0
    as.matrix(cm) = cm_matrix
    iset = deflate(cm, use.zero = TRUE)
    iset
    })
  dense_iset = do.call(cbind, dense_isets)
  assayNames(dense_iset) = 'counts'
  dense_iset
}

filter_regions_iset <- function(iset, filter_regions_path, resolution){  
  filter_regions = fread(filter_regions_path)
  filter_regions = filter_regions[,c(1:2)]
  colnames(filter_regions) = c('chr','start')
  filter_regions$bin_id = as.integer(filter_regions$start/resolution)
  #filter_regions_dt = as.data.table(filter_regions)
  #filter_regions_dt = unique(filter_regions_dt, by = c('chr','bin_id'))
  #filter_regions_dt$include = rep(FALSE, length(filter_regions_dt$chr))
  
  filter_regions = data.frame(filter_regions %>% group_by(chr,bin_id) %>% summarise(n = n()))
  
  filter_regions = filter_regions[filter_regions$n == as.integer(resolution/10000),]
  #filter_regions = filter_regions[filter_regions$n > 0,]
  filter_regions_dt = as.data.table(filter_regions)
  filter_regions_dt$include = rep(FALSE, length(filter_regions_dt$chr))
  
  binIds <- anchors(iset, type="both")
  binId_dt <- data.table(bin1_id = as.integer((data.frame(binIds$first)$start)/resolution), 
                         bin2_id= as.integer((data.frame(binIds$second)$start)/resolution), 
                         chrs= data.frame(binIds$first)$seqnames) 
  # Perform the join and filtering operation
  binId_dt = left_join(binId_dt, filter_regions_dt, by = c("bin1_id" = "bin_id", "chrs" ="chr"))
  binId_dt = left_join(binId_dt, filter_regions_dt, by = c("bin2_id" = "bin_id", "chrs" ="chr"))
  binId_dt[is.na(include.x), include.x := TRUE]
  binId_dt[is.na(include.y), include.y := TRUE]
  
  keep_filter_regions = binId_dt$include.x & binId_dt$include.y
  iset <- iset[keep_filter_regions,]
  
  return(iset)
}
