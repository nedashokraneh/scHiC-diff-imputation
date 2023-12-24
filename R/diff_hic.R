library(diffHic)
library(csaw)
library(edgeR)


get_diffhic_result_table <- function(iset, group){
  y <- asDGEList(iset)
  y <- normOffsets(y, se.out=TRUE)
  design = model.matrix(~0+group)
  contrast12 <- makeContrasts( Group2vs1 = group2-group1,
                               levels = design )
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust=TRUE)
  diffhic_result <- glmQLFTest(fit, contrast = contrast12)$table
  rg1 = data.frame(anchors(iset)$first)$start 
  rg2 = data.frame(anchors(iset)$second)$start 
  chroms = data.frame(anchors(iset)$first)$seqnames
  
  modify_res <- function(res){
    names(res)[names(res) == 'PValue'] <- 'p.value'
    res$p.adj <- p.adjust(res$p.value, method="BH")
    res$chrom = chroms
    res$region1 <- rg1
    res$region2 <- rg2
    res
  }
  diffhic_result = modify_res(diffhic_result)
  return(diffhic_result)
}