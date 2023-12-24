library('rhdf5')
library(pryr)
#library(purrr)
library(InteractionSet)
library(dplyr)
library(data.table)
library(Matrix)
library(sparsesvd)
library(umap)
library(csaw)
library(edgeR)
library(multiHiCcompare)
library(pbapply)
#TODO: add txt loader

### this function is from https://github.com/TaoYang-dev/hicrep/blob/master/R/cool2matrix.R
cool2matrix <- function(file, chr = 'chr1') {

  pixels <- h5read(file, c('pixels'));H5close()
  bins <- h5read(file, c('bins'));H5close()
  chrom.offset <- h5read(file, 'indexes/chrom_offset');H5close()

  chr2index <- function(chr) {
    if (substr(chr, 1, 3) == 'chr') {
      chr = substr(chr, 4, nchar(chr))
    }
    if (chr == 'X') {
      index = 23
    }
    else if (chr == 'Y') {
      index = 24
    }
    else if (chr == 'M') {
      index = 25
    }
    else {
      index = as.integer(chr)
    }
    return(index)
  }

  chrom.index <- chr2index(chr)
  chrom.range <- chrom.offset[chrom.index:(chrom.index+1)] + c(0, -1) 
  n.rows <- chrom.range[2] - chrom.range[1] + 1 
  bin.range <- which(pixels$bin1_id >= chrom.range[1] & 
                       pixels$bin1_id <= chrom.range[2] & 
                       pixels$bin2_id >= chrom.range[1] & 
                       pixels$bin2_id <= chrom.range[2])
  n.bins <- length(bin.range) 
  mat <- matrix(0, ncol = n.rows, nrow = n.rows)
  for (i in 1:n.bins) {
    mat[pixels$bin1_id[bin.range[i]] - chrom.range[1] + 1, 
        pixels$bin2_id[bin.range[i]] - chrom.range[1] + 1] <- 
      pixels$count[bin.range[i]]
    mat[pixels$bin2_id[bin.range[i]] - chrom.range[1] + 1, 
        pixels$bin1_id[bin.range[i]] - chrom.range[1] + 1] <- 
      pixels$count[bin.range[i]]
  }
  chrom.ranges <- (chrom.range[1]+1):(chrom.range[2]+1)
  # mat <- as.data.frame(mat)
  
  return(mat)
}


cool2sparse <- function(file, chr = 'chr1', resolution) {
  pixels <- h5read(file, c('pixels'));H5close()
  bins <- h5read(file, c('bins'));H5close()
  chrom.offset <- h5read(file, 'indexes/chrom_offset');H5close()
  chr2index <- function(chr) {
    if (substr(chr, 1, 3) == 'chr') {
      chr = substr(chr, 4, nchar(chr))
    }
    if (chr == 'X') {
      index = 23
    }
    else if (chr == 'Y') {
      index = 24
    }
    else if (chr == 'M') {
      index = 25
    }
    else {
      index = as.integer(chr)
    }
    return(index)
  }

  chrom.index <- chr2index(chr)
  chrom.range <- chrom.offset[chrom.index:(chrom.index+1)] + c(0, -1) 
  df = data.frame(bin1_id = pixels$bin1_id,
                  bin2_id = pixels$bin2_id, 
                  count = pixels$count)
  df = df[(df$bin1_id >= chrom.range[1]) & 
            (df$bin1_id <= chrom.range[2]) & 
            (df$bin2_id >= chrom.range[1]) & 
            (df$bin2_id <= chrom.range[2]),]
  df[,c('bin1_id','bin2_id')] = (df[,c('bin1_id','bin2_id')]-chrom.range[1])*resolution
  df$chr = chr
  df = df[,c('chr','bin1_id','bin2_id','count')]
  return(df)
}

### convert 3 columns data frame to contact map
df2mat <- function(df, chr_size, resolution){
  colnames(df) = c('bin1_id', 'bin2_id', 'count')
  chr_size = ceiling(chr_size/resolution)
  mat <- matrix(0, nrow=chr_size, ncol=chr_size)
  df[,c('bin1_id','bin2_id')] = df[,c('bin1_id','bin2_id')] / resolution
  mat[as.matrix(df[,c('bin1_id','bin2_id')]+1)] <- df[,'count']
  mat[as.matrix(df[,c('bin2_id','bin1_id')]+1)] <- df[,'count']
  mat
}

make_diffhic_object <- function(input_format, files, chr_name, 
                                chr_size_file, resolution){
  g = read.table(chr_size_file)
  chr_size = as.numeric(g[g[,1]==chr_name,2])
  if (input_format == 'txt'){
    dfs = lapply(files, read.table)
    mats = lapply(dfs, df2mat, chr_size = chr_size, resolution = resolution)
  }
  else if (input_format == 'cool'){
    mats = lapply(files, cool2matrix, chr = chr_name)
  }
  else{
    print('txt or cool formats are available.')
    return 
  }
  print(object_size(mats))
  chr_size = ceiling(chr_size/resolution)
  regions = GRanges(rep(chr_name,chr_size), IRanges(c(0:(chr_size-1)),c(1:chr_size)))
  cms = lapply(mats, ContactMatrix, anchor1 = c(1:chr_size),
                     anchor2 = c(1:chr_size), regions = regions)
  print(object_size(cms))
  to.keep = Reduce("|", lapply(cms, function(cm){as.matrix(cm)!=0}))
  isets = lapply(cms, deflate, extract = to.keep)
  data = Reduce(cbind, isets)
  interactions(data) <- as(interactions(data), "ReverseStrictGInteractions")
  assayNames(data) = 'counts'
  return(data)
}

# make_multiHiCcompare_object <- function(input_format, files, chr_name, 
#                                         chr_size_file, resolution, groups){
#   if (input_format == 'txt'){
#     dfs = lapply(files, function(f){
#       df = read.table(f)
#       colnames(df) = c('bin1_id', 'bin2_id', 'count')
#       df$chr = chr_name
#       df = df[,c('chr', 'bin1_id', 'bin2_id', 'count')]
#       return (df)
#     })
#   }
#   else if (input_format == 'cool'){
#     dfs = lapply(files, cool2sparse, chr = chr_name, resolution = resolution)
#   }
#   else{
#     print('txt or cool formats are available.')
#     return 
#   }
#   hicexp = make_hicexp(data_list=dfs, groups = groups, A.min = 5, zero.p = 1)
#   return(hicexp)
# }

bin_hic <- function(df, resolution){
  df$bin1_id = as.integer(df$bin1_id/resolution)
  df$bin2_id = as.integer(df$bin2_id/resolution)
  min_bin_id = pmin(df$bin1_id, df$bin2_id)
  max_bin_id = pmax(df$bin1_id, df$bin2_id)
  df$bin1_id = min_bin_id
  df$bin2_id = max_bin_id
  d_threshold = max(1000000/resolution,10)
  df = df[df$bin2_id-df$bin1_id <= d_threshold,]
  binned_df <- setDT(df)[,list(count=.N),names(df)]
  binned_df
}

make_regions_vecs <- function(regions_df, include, chrom_size_filepath, 
                                chromosomes, resolution){
  chrom_sizes = read.table(chrom_size_filepath)
  bool_vecs = list()
  for (chrom in chromosomes){
    region_ids = regions_df[regions_df$chr == chrom, 'bin_id']
    region_ids = as(region_ids, 'vector')$bin_id
    chrom_size = ceiling(chrom_sizes[chrom_sizes$V1 == chrom, 'V2']/resolution)
    bool_vecs[[chrom]] = vector(mode = "logical", length = chrom_size)
    bool_vecs[[chrom]][1:chrom_size] = !include
    bool_vecs[[chrom]][region_ids] = include
  }
  bool_vecs
}

### TODO: faster filtering? filtering before or after normalization?

read_hic_df <- function(path, chrom_columns, pos_columns, resolution,
                        filter_regions_path, TSS_regions_path, chrom_size_filepath,
                        chromosomes){
  df = data.frame(fread(path, sep = "\t"))
  df = df[df[,chrom_columns[1]]==df[,chrom_columns[2]],]
  #df = df[df[,chrom_columns[1]] %in% paste0('chr', c(1:22)),]
  df = df[df[,chrom_columns[1]] %in% chromosomes,]
  df = df[,c(chrom_columns[1],pos_columns)]
  colnames(df) = c('chr', 'bin1_id', 'bin2_id')
  binned_df = bin_hic(df,resolution)
  
  # only including valid bins 
  # filter_regions = fread(filter_regions_path)
  # filter_regions = filter_regions[,c(1:2)]
  # colnames(filter_regions) = c('chr','start')
  # filter_regions$bin_id = as.integer(filter_regions$start/resolution)
  # filter_vecs = make_regions_vecs(filter_regions, FALSE, chrom_size_filepath, 
  #                                 paste0('chr', c(1:22)), resolution)
  # is_valid = function(x){filter_vecs[[x['chr']]][as.integer(x['bin1_id'])] & 
  #     filter_vecs[[x['chr']]][as.integer(x['bin2_id'])]}
  # valid = apply(binned_df, 1, is_valid)
  # valid = as.numeric(valid)
  # valid[is.na(valid)] = 0
  # 
  # # only including TSS regions 
  # TSS_regions = fread(TSS_regions_path)
  # TSS_regions$bin_id = as.integer((TSS_regions$start+TSS_regions$end)/(2*resolution))
  # TSS_vecs = make_regions_vecs(TSS_regions, TRUE, chrom_size_filepath, 
  #                                 paste0('chr', c(1:22)), resolution)
  # #return(list('TSS_vecs'=TSS_vecs,'df'=binned_df))
  # is_TSS = function(x){TSS_vecs[[x['chr']]][as.integer(x['bin1_id'])] |
  #     TSS_vecs[[x['chr']]][as.integer(x['bin2_id'])]}
  # TSS = apply(binned_df, 1, is_TSS)
  # TSS = as.numeric(TSS)
  # TSS[is.na(TSS)] = 0
  # binned_df = binned_df[valid&TSS,]
  binned_df
}

read_hic_cool <- function(cond1_cooldir, cond2_cooldir, chr, sample_size){
  coolpaths = lapply(list(cond1_cooldir, cond2_cooldir), function(p){
    file.path(p,list.files(p))})
  coolpaths[[1]] = sample(coolpaths[[1]], sample_size)
  coolpaths[[2]] = sample(coolpaths[[2]], sample_size)
  groups = unlist(lapply(c(1:2), function(i){
    rep(paste0('cond', i), length(coolpaths[[i]]))}))
  coolpaths = unlist(coolpaths)
  all_pixels = pblapply(c(1:length(coolpaths)), function(coolpaths, id){print(id)
    pixels <- h5read(coolpaths[id], c('pixels'));H5close()
    pixels = data.frame(pixels)
    pixels$chr = chr
    pixels[c('chr', 'bin1_id', 'bin2_id', 'count')]
  }, coolpaths = coolpaths)
  all_pixels
}

# input: list of hic dataframes including chr, bin1_id, bin2_id, and count columns
# output: one hic dataframe including all cells 

rbind_hic_dfs <- function(hic_dfs){
  all_hic_df = lapply(c(1:length(hic_dfs)), function(id){
    df = hic_dfs[[id]]
    df = df[,c('chr','bin1_id','bin2_id','count')]
    df$cell = paste0('cell',id)
    df
  })
  all_hic_df = rbindlist(all_hic_df)
  all_hic_df
}

# input: list of hic dataframes including chr, bin1_id, bin2_id, and count columns
# input: resolution
# input: vector of cell types corresponding to cells
# output: hicexp object for multiHiCcompare
# TODO: zero.p is used for single cell resolution, A.min for pseudobulk?

make_hicexp_obj <- function(hic_dfs, res, groups){
  hicexp_dfs = lapply(hic_dfs, function(df){
    df = df[df$bin1_id!=df$bin2_id,]
    df[,c(2,3)] = df[,c(2,3)]*res
    data.frame(df)
  })
  hicexp = make_hicexp(data_list=hicexp_dfs, groups = groups, A.min = 0, zero.p = 1)
}

# input: one hic dataframe including chr, bin1_id, bin2_id, count, and cell columns
# output: one hic dataframe with a shape of original one including BandNorm counts

bandnorm <- function(hic_df){
  hic_df$diag = hic_df$bin2_id - hic_df$bin1_id
  hic_df = hic_df[hic_df$diag != 0, ]
  band_info <- hic_df %>% group_by(diag, cell) %>% summarise(band_depth = sum(count))
  alpha_j <- band_info %>% group_by(diag) %>% summarise(depth = mean(band_depth))
  hic_df <- hic_df %>% left_join(alpha_j, by = c("diag")) %>% 
    left_join(band_info, by = c("diag", "cell")) %>% 
    mutate(BandNorm = count/band_depth * depth) %>% 
    select(-c(band_depth, depth, count))
  #colnames(hic_df) = c('bin1_id', 'bin2_id', 'cell', 'diag', 'count')
  hic_df
}

# input: list of hic dataframes including chr, bin1_id, bin2_id, and count columns
# output: one InteractionSet object 

make_iset <- function(dfs){
  dfs = lapply(c(1:length(dfs)), function(data_frames, id){
    df = data_frames[[id]]
    colnames(df) = c('chr', 'bin1_id', 'bin2_id', paste0('count',id))
    data.table(df, key = c('chr', 'bin1_id', 'bin2_id'))
  }, data_frames = dfs)
  merged_df =  Reduce(function(dt1,dt2){merge.data.table(dt1,dt2,all = TRUE)}, 
                      dfs)
  merged_df[is.na(merged_df)] = 0
  chrs = merged_df$chr
  source_gr = GRanges(seqnames = merged_df$chr, 
                      ranges = IRanges(merged_df$bin1_id, end = merged_df$bin1_id+1))
  target_gr = GRanges(seqnames = merged_df$chr,
                      ranges = IRanges(merged_df$bin2_id, end = merged_df$bin2_id+1))
  GI = GInteractions(source_gr, target_gr)
  iset = InteractionSet(as(merged_df[,4:ncol(merged_df)], 'matrix'), GI)
  elementMetadata(iset) = chrs
  iset
}

### TODO: coarsen_iset?

# input: one hic dataframe including all cells with chr, bin1_id, bin2_id, \
# count, and cell columns or one hic matrix of size #cells \times #IFs
# input: format is df or mat
# input: vector of cell types corresponding to cells
# output: dataframe with a cell per row and columns correponding to its UMAP embeddings and cell type 

hic_embedding <- function(hic_data, format, groups){
  if (format == 'mat'){
    hic_mat = as(hic_data, 'sparseMatrix')
    #hic_mat = hic_mat[,colSums(hic_mat!=0)>20]
  }
  if (format == 'df'){
    hic_df = hic_data
    hic_df$diag = hic_df$bin2_id - hic_df$bin1_id
    hic_df = hic_df[hic_df$diag != 0, ]
    names = unique(hic_df$cell)
    hic_df$cell = as.numeric(factor(hic_df$cell, level = names))
    hic_df = hic_df %>% 
      mutate(featureIndex = paste(bin1_id, bin2_id, sep = "_")) %>%
      select(-c(bin1_id, bin2_id, diag))
    hic_df$featureIndex = as.numeric(factor(hic_df$featureIndex))
    hic_mat = sparseMatrix(i = hic_df$cell, j = hic_df$featureIndex, x = hic_df$count, 
                           dims = c(max(hic_df$cell), max(hic_df$featureIndex)),
                           index1 = TRUE)
    hic_mat = hic_mat[,colSums(hic_mat!=0)>dim(hic_mat)[2]/10]
  }
  
  pca_mat = sparsesvd(hic_mat, 20)
  pca_mat = pca_mat$u %*% diag(pca_mat$d)
  
  embedding = umap(pca_mat)$layout
  embedding = data.frame(embedding)
  colnames(embedding) = c("X1", "X2")
  embedding$cell_type = groups
  embedding
}



# input: one InteractionSet object
# output: one InteractionSet that has been filtered to 
# exclude filter regions, and only include TSS regions


filter_regions_df <- function(df, filter_regions_path, resolution){  
  filter_regions = fread(filter_regions_path)
  filter_regions = filter_regions[,c(1:2)]
  colnames(filter_regions) = c('chr','start')
  filter_regions$bin_id = as.integer(filter_regions$start/resolution)
  filter_regions_dt = as.data.table(filter_regions)
  filter_regions_dt = unique(filter_regions_dt, by = c('chr','bin_id'))
  filter_regions_dt$include = rep(FALSE, length(filter_regions_dt$chr))
  
  binId_dt = data.table(bin1_id = df$bin1_id, 
                        bin2_id = df$bin2_id, 
                        chr = df$chr) 
  # Perform the join and filtering operation
  binId_dt = left_join(binId_dt, filter_regions_dt, by = c("bin1_id" = "bin_id", "chr" ="chr"))
  binId_dt = left_join(binId_dt, filter_regions_dt, by = c("bin2_id" = "bin_id", "chr" ="chr"))
  binId_dt[is.na(include.x), include.x := TRUE]
  binId_dt[is.na(include.y), include.y := TRUE]
  
  keep_filter_regions = binId_dt$include.x & binId_dt$include.y
  df <- df[keep_filter_regions,]
  
  return(df)
}


# input: one InteractionSet object
# output: one InteractionSet that has been filtered to 
# exclude filter regions, and only include TSS regions
TSS_filter_iset <- function(iset, TSS_regions_path, resolution){  
  # # only including TSS regions 
  TSS_regions = fread(TSS_regions_path)
  TSS_regions$bin_id = as.integer((TSS_regions$start+TSS_regions$end)/(2*resolution))
  TSS_regions_dt = as.data.table(TSS_regions)
  TSS_regions_dt$include = rep(TRUE, length(TSS_regions_dt$chr))
  TSS_regions_dt = unique(TSS_regions_dt[,c("bin_id","include","chr")])
  
  chrs <- elementMetadata(iset)$X
  binIds <- anchors(iset, type="both")
  binId_dt <- data.table(bin1_id = data.frame(binIds$first)$start, 
                         bin2_id= data.frame(binIds$second)$start, chrs=chrs) 
  
  # Perform the join and filtering operation
  binId_dt = left_join(binId_dt, TSS_regions_dt, by = c("bin1_id" = "bin_id", "chrs" ="chr"))
  binId_dt = left_join(binId_dt, TSS_regions_dt, by = c("bin2_id" = "bin_id", "chrs" ="chr"))
  binId_dt[is.na(include.x), include.x := FALSE]
  binId_dt[is.na(include.y), include.y := FALSE]
  
  keep_gene_transcript = binId_dt$include.x & binId_dt$include.y
  iset <- iset[keep_gene_transcript,]

  return(iset)
}

diag_filter_iset <- function(iset){
  binIds <- anchors(iset, type="both")
  keep = data.frame(binIds$first)$start != data.frame(binIds$second)$start
  iset = iset[keep,]
  iset
}

### TODO: coarsen_iset?
# input: one InteractionSet
# input: n_groups, the number pseudo_bulk groups per cell type
coarsen_iset <- function(isets, n_groups){
  n_samples = dim(assay(isets))[2]
  if ((n_samples/2) %% n_groups != 0){
    print("pseudo-bulk will not have the same number of samples")
    return(NULL)
  }
  
  group_size = (n_samples/2) / n_groups
  for (i in seq(1, n_samples - group_size + 1, by = group_size)) {
    new_col <- rowSums(assay(isets)[, i:(i + group_size - 1)])
    col_name <- paste0("coarse_", i, "-", i + group_size - 1)
    col_names = colnames(isets)
    col_names[ceiling(i/group_size)]=col_name
    colnames(isets)=col_names
    assay(isets)[, ceiling(i/group_size)]= new_col
  }
  
  # Remove the original columns
  isets = isets[, 1:(n_samples/group_size)]
  return(isets)
}

calc_iset_sparsity <- function(iset){
  cnts = assays(iset)[[1]]
  nz_cnts = sum(cnts != 0)
  nz_cnts / (dim(cnts)[1] * dim(cnts)[2])
}

prep_isets <- function(dfs, filter_regions_path, experiment_sizes){
  iset = make_iset(dfs)
  iset = diag_filter_iset(iset)
  max_es = experiment_sizes[length(experiment_sizes)]
  iset = iset[,c(1:max_es, (dim(iset)[2]-max_es+1):dim(iset)[2])]
  assayNames(iset) = 'counts'
  iset = filter_regions_iset(iset, filter_regions_path, 100000)
  keep <- (rowMeans(assay(iset)[,1:max_es]) > 0.1) | 
    (rowMeans(assay(iset)[,(max_es+1):(2*max_es)]) > 0.1)
  iset = iset[keep,]
  isets = list()
  for (es in experiment_sizes[1:(length(experiment_sizes)-1)]){
    isets[[es]] = coarsen_iset(iset, es)
  }
  isets[[max_es]] = iset
  isets
}

