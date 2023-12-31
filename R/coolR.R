# adapted from https://github.com/dovetail-genomics/coolR/blob/master/R/coolR.R

#' Accessing the sequence information in a cool/mcool file
#'
#' Convering the chromosomes information from a HDF5 stored uni- or multi-dimentional
#' contact matrix (cool or mcool file) into a seqinfo object
#' 
#' @param file Path to a HDF5 stored cool (uni-dimension) or mcool (multi-dimension sparse matrix).
#' If using a multi-dimesion, the resolution of one of the dimension need to be passed to res
#' @param res 'NULL' if using a uni-dimensional cool file or the resolution of one of layer in the mcool file
#'
#' @examples
#' 
#' @export
seqinfo.cool <- function(file,res = NULL){
  
  chr.group <- ifelse(is.null(res),'/chroms',paste('resolutions',res,'chroms',sep='/'))
  
  s.i <- Seqinfo(seqnames = as.vector(h5read(file,name=paste(chr.group,'name',sep="/"))),
                 seqlengths = as.vector(h5read(file,name=paste(chr.group,'length',sep="/"))))
  
  return(s.i)
}

#' Accessing the genomic bins in a cool/mcool file
#'
#' Retrieves the chromosome bins forming from a HDF5 stored uni- or multi-dimentional
#' contact matrix (cool or mcool file) into a seqinfo object
#' 
#' @param file Path to a HDF5 stored cool (uni-dimension) or mcool (multi-dimension sparse matrix).
#' If using a multi-dimesion, the resolution of one of the dimension need to be passed to res
#' @param res 'NULL' if using a uni-dimensional cool file or the resolution of one of layer in the mcool file
#'
#' @examples
#' 
#' @export
getBins <- function(file,res = NULL){
  
  bins.group <- ifelse(is.null(res),'/bins',paste('resolutions',res,'bins',sep='/'))
  
  s.i <- seqinfo.cool(file,res)
  
  a.chr <- h5read(file,name=paste(bins.group,'chrom',sep="/"))
  ## Cooler use open ended starts, we will use close non-overlapping bins
  a.start <- h5read(file,name=paste(bins.group,'start',sep="/"))+1
  a.end <- h5read(file,name=paste(bins.group,'end',sep="/"))
  
  anchors <- GRanges(a.chr,
                     IRanges(a.start,a.end),
                     seqinfo=s.i)
  
  return(anchors)
}

getSlice <- function(anchors,file,res,chr1,start1,end1,chr2,start2,end2){
  ################################################################################
  ## Get the chromosme bins as a GRanges
  ################################################################################
  if(is.null(start1)) start1 <- 1
  if(is.null(end1)) end1 <- seqlengths(anchors)[chr1]
  
  indexes.group <- ifelse(is.null(res),"/indexes",paste('resolutions',res,'indexes',sep='/'))
  chr.group <- ifelse(is.null(res),"/chroms",paste('resolutions',res,'chroms',sep='/'))
  pixels.group <- ifelse(is.null(res),"/pixels",paste('resolutions',res,'pixels',sep='/'))
  
  indexes <- list(chr=as.vector(h5read(file,name=paste(chr.group,'name',sep="/"))),
                  chr_idx=as.vector(h5read(file,name=paste(indexes.group,'chrom_offset',sep="/"))),
                  bin1_idx=as.vector(h5read(file,name=paste(indexes.group,'bin1_offset',sep="/")))
  )
  
  ################################################################################
  ## Reading chromosome chunk If chr1 is null, return the full cool file
  ################################################################################
  if (is.null(chr1)){
    chunk <- NULL
  } else {
    chr1.start <- GRanges(chr1,IRanges(start1,width=1))
    chr1.end <- GRanges(chr1,IRanges(end1,width=1))
    
    chr1.start.idx <- subjectHits(findOverlaps(chr1.start,anchors))
    chr1.end.idx <- subjectHits(findOverlaps(chr1.end,anchors))
    idx.chunk <- seq(chr1.start.idx,chr1.end.idx)
    
    
    bin1.idx <- as.vector(h5read(file,
                                 name=paste(indexes.group,'bin1_offset',sep="/"),
                                 index=list(idx.chunk)))
    
    slice <- sum(bin1.idx[-1] - bin1.idx[-length(bin1.idx)])-1
    
    chunk  <- seq(bin1.idx[1]+1,bin1.idx[1]+1+slice)
    
  }
  
  ################################################################################
  ## Reading the chunks from the cool file
  ################################################################################
  d.f <- data.frame(bin1_id = as.vector(h5read(file,
                                               name=paste(pixels.group,'bin1_id',sep="/"),
                                               index=list(chunk)))+1,
                    bin2_id = as.vector(h5read(file,
                                               name=paste(pixels.group,'bin2_id',sep="/"),
                                               index=list(chunk)))+1,
                    count = as.vector(h5read(file,
                                             name=paste(pixels.group,'count',sep="/"),
                                             index=list(chunk)))
  )
  
  ################################################################################
  ## If Chr2 is set, only return corresponding ranges
  ################################################################################
  if (!is.null(chr2)){
    if(is.null(start2)) start2 <- 1
    if(is.null(end2)) end2 <- seqlengths(anchors)[chr2]
    
    chr2.start <- GRanges(chr2,IRanges(start2,width=1))
    chr2.end <- GRanges(chr2,IRanges(end2,width=1))
    chr2.start.idx <- subjectHits(findOverlaps(chr2.start,anchors))
    chr2.end.idx <- subjectHits(findOverlaps(chr2.end,anchors))
    filter.bin2 <- d.f$bin2_id %in% chr2.start.idx:chr2.end.idx
    
    d.f <- d.f[filter.bin2,]
  }
  
  chr1_offset = indexes$chr_idx[which(indexes$chr == chr1)]
  chr2_offset = indexes$chr_idx[which(indexes$chr == chr2)]
  d.f$bin1_id = d.f$bin1_id - chr1_offset
  d.f$bin2_id = d.f$bin2_id - chr2_offset
  
  return(d.f)
  
}


gi2is <- function(gi.counts,col.names) {
  gi <- unique(do.call(c,gi.counts))    
  S4Vectors::mcols(gi) <- c() 
  
  ## Get a dataframe of counts for each interaction pairs
  counts <- do.call(cbind,lapply(gi.counts,function(x){
    ## Initialize a vector of counts 0 of the lenght of all interaction
    counts <- as.integer(rep(0,length(gi)))
    ## Find the interactions in the common set intersecting with the current GI
    mm <- findMatches(gi,x)
    ## Replace in the vector of counts for the universe, replace the intersection with the current counts
    counts[queryHits(mm)] <- x$count[subjectHits(mm)]
    ## Return
    return(counts)
  }))
  
  ## Asign the name of the Column to the sample
  colnames(counts) <- col.names
  
  ## Get the library size
  lib.data <- DataFrame(totals=colSums(counts))
  
  ## Create the interaction set
  data <- InteractionSet(assays=list(counts=counts),
                         gi,
                         colData=lib.data)
  
  return(data)
}


#' Read .cool/.mcool sparse matrix into a GInteractions Object
#'
#' This function reads a compressed sparse row (CSR) storage scheme for a matrix
#' created by the cooler application (https://github.com/mirnylab/cooler) stored in a
#' HDF5 data storage. This utility can read the matrix or only specfic region of the genome
#' eihter from one pair or a pair of region. If using a multi-dimension dataset (ie .mcool), the resolution
#' needed need to be specified.
#'
#' @param file Path to a HDF5 stored cool (uni-dimension) or mcool (multi-dimension sparse matrix).
#' If using a multi-dimesion, the resolution of one of the dimension need to be passed to res
#' @param res 'NULL' if using a uni-dimensional cool file or the resolution of one of layer in the mcool file
#' @param chr1 'NULL' or the chromosome name of the first pairs to return
#' @param start1 'NULL' or the start of a region within chr1 to return the pairs from.
#' If chr1 is set and start1 is NULL this will be set to 1 by default
#' @param end1 'NULL' or the end of a region within chr1 to return the pairs from.
#' If chr1 is set and end1 is NULL this will be set to the lenght of chr1  by default
#' @param chr2 'NULL' or the chromosome name of the first pairs to return
#' @param start2 'NULL' or the start of a region within chr2 to return the pairs from.
#' If chr1 is set and start1 is NULL this will be set to 1 by default
#' @param end2 'NULL' or the end of a region within chr1 to return the pairs from.
#' If chr2 is set and end1 is NULL this will be set to the lenght of chr2  by default
#'
#' @return A GInteraction object of the contact matrix dataset
#' 
#' @import InteractionSet
#' @import rhdf5
#' 
#' @export
read.cool <- function(file,res=NULL,chr1=NULL,start1=NULL,end1=NULL,chr2=NULL,start2=NULL,end2=NULL) {
  
  anchors <- getBins(file,res)
  
  slice <- getSlice(anchors,file,res,chr1,start1,end1,chr2,start2,end2)
  
  # g.i <- GInteractions(anchors[slice$bin1_id],
  #                      anchors[slice$bin2_id],
  #                      count=slice$count)
  
  all.regions = anchors[seqnames(anchors) %in% c(chr1, chr2)]
  
  if (chr1 != chr2){
    chr1_num_regions = length(anchors[seqnames(anchors) == chr1])
    slice$bin2_id = slice$bin2_id + chr1_num_regions
  } 

  g.i <- GInteractions(slice$bin1_id,
                       slice$bin2_id,
                       all.regions,
                       count=slice$count)
  
  return(g.i)
}


#' Read .cool/.mcool sparse matrix into a InteractionSet object
#'
#' This function reads HiC contact matrix file(s) created by the 
#' cooler application (https://github.com/mirnylab/cooler) stored in a
#' HDF5 data storage and converts them to an InteractionSet object for all
#' detected pair of genomic bins with the counts of reads in these bin pairs.
#' Reading and processing the files are done in parallel using mclapply from the
#' parallel package and usually requires ~1Gb of RAM per core.
#' 
#' @param file.path The path to a HDF5 stored cool (uni-dimension)
#' or mcool (multi-dimension sparse matrix). If using a multi-dimesion,
#' the resolution of one of the dimension need to be passed to res.
#' @param ... additional path file(s) to HDF5 stored cool (uni-dimension)
#' or mcool (multi-dimension sparse matrix).
#' @param res 'NULL' if using a uni-dimensional cool file or the resolution of one of layer in the mcool file
#' @param cores An integer for the number of parallel thread to convert the file
#'
#' @return An InteractionSet of all common pairs with the number of interaction for each pairs
#' 
#' @import InteractionSet
#' @import rhdf5
#'
#' @export
cool2IntSet <- function(file.path, ..., res=NULL, cores = detectCores() ) {
  
  files <- c(file.path,...)
  
  anchors <- lapply(files,getBins,res)
  
  if (!(all(sapply(lapply(anchors[-1],seqinfo),identical,seqinfo(anchors[[1]]))))){
    stop('cooler files do not have same chromosomes')
  }
  
  map <- data.frame(cool.f = rep(files,sapply(anchors,function(x) length(seqlevels(x)))),
                    chr = do.call(c,lapply(anchors,seqlevels)))
  
  gi.l <- mclapply(seq(nrow(map)),function(i){
    chr <- map$chr[i]
    file <- map$cool.f[i]
    read.cool(file,res,chr1=chr)
  },mc.cores=cores,mc.preschedule=FALSE)
  
  
  ## Run the reduce step
  suppressWarnings({
    res <- do.call(c, mclapply(unique(map$chr),function(chr){
      ints <- gi2is(gi.l[map$chr == chr],
                    col.names = basename(files))
    },mc.cores=cores,mc.preschedule=FALSE))
  })
  
  return(res)
}

