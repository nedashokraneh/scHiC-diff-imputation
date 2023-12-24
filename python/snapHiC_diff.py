# this script is based on https://github.com/HuMingLab/SnapHiC-D/blob/main/SnapHiC_D/SnapHiC_D.py

import argparse

import qnorm
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats
import statsmodels.stats.multitest as multi


tstat_th = 2
fdr_th = 0.05

def main():
    
    parser = create_parser()
    args = parser.parse_args()
    
    if args.batch is None:
        comp_batch_names = None
    else:
        batch1 = args.batch
        if args.batch2 is None:
            batch2 = batch1
            
        else:
            batch2 = args.batch2 
        comp_batch_names = [batch1, batch2]
            
    if 'downsample' in args:
        downsample = args.downsample
    else:
        downsample = None 
        
    comp_cell_types = [args.groupA, args.groupB]
    
    all_anndata_paths = read_anndata_filelist(args.filelist)
    empty_out = True
    for anndata_path in all_anndata_paths.keys():
        chromosomes = all_anndata_paths[anndata_path]
        adata = ad.read_h5ad(anndata_path)
        if comp_batch_names is None:
            #keep_comp_cells = [ct in comp_cell_types for ct in adata.obs["cell_type"]]
            idx1 = np.where(adata.obs["cell_type"] == comp_cell_types[0])[0]
            idx2 = np.where(adata.obs["cell_type"] == comp_cell_types[1])[0]
        
        else:
            #keep_comp_cells = [((ct == comp_cell_types[0]) and (b == batch1)) or ((ct == comp_cell_types[1]) and (b == batch2)) for ct, b in zip(adata.obs["cell_type"], adata.obs["batch_name"])]
            idx1 = np.where((adata.obs["cell_type"] == comp_cell_types[0]) & (adata.obs["batch_name"] == comp_batch_names[0]))[0]
            idx2 = np.where((adata.obs["cell_type"] == comp_cell_types[1]) & (adata.obs["batch_name"] == comp_batch_names[1]))[0]
        
        if downsample is not None:
            idx1 = np.random.choice(idx1, downsample, replace = False)
            idx2 = np.random.choice(idx2, downsample, replace = False)
        
        idx = np.concatenate([idx1, idx2])
        adata = adata[idx, :]
        keep_nz_feats = adata.X.sum(axis = 0) != 0
        adata = adata[:, keep_nz_feats]
        #keep_var_feats = np.var(adata.X.toarray(), ddof = 1, axis = 0) != 0
        #adata = adata[:, keep_var_feats]
        
        for chromosome in chromosomes:
            adata_ = adata[:, adata.var["chrom"] == chromosome].copy()
            keep_var_feats = np.var(adata_.X.toarray(), ddof = 1, axis = 0) != 0
            adata_ = adata_[:, keep_var_feats]
            adata_, norm_layer = normalize(adata_, args.norm)
            mask1, mask2 = get_groups_masks(adata_, comp_cell_types, comp_batch_names)
            if norm_layer == 'X':
                A_scores = adata_.X[mask1, :].toarray()
                B_scores = adata_.X[mask2, :].toarray()
            else:
                A_scores = adata_.layers[norm_layer][mask1, :] #.toarray()
                B_scores = adata_.layers[norm_layer][mask2, :] #.toarray()
            
            
            t_stat, p_val = stats.ttest_ind(np.log2(A_scores+1), np.log2(B_scores+1), axis=0, equal_var=False)
            real_p_val_idx = np.where(~np.isnan(p_val))[0]
            p_val_ = p_val[real_p_val_idx]
            fdr_ = multi.multipletests(p_val_, method='fdr_bh')[1]
            fdr = np.ones(p_val.shape)
            fdr[real_p_val_idx] = fdr_
            result_table = pd.DataFrame(adata_.var).copy()
            result_table['mean.A'] = np.array(np.mean(A_scores, axis=0)) 
            result_table['mean.B'] = np.array(np.mean(B_scores, axis=0))
            result_table['Tstat'] = t_stat
            result_table['Ttest.Pvalue'] = p_val
            result_table['fdr'] = fdr
            result_table['significant'] = (result_table['fdr'] < fdr_th) & (abs(result_table['Tstat']) > tstat_th)
            if empty_out:
                result_table.to_csv(args.output_path, sep = "\t", index = False, mode = 'w')
                empty_out = False
            else:
                result_table.to_csv(args.output_path, sep = "\t", header = False, index = False, mode = 'a')

            
'''
def get_groups_masks(adata, cell_types, batch_name = None):
    assert len(cell_types) == 2, "two cell types need to be specified."
    if batch_name is not None:
        mask1 = (adata.obs["cell_type"] == cell_types[0]) & (adata.obs["batch_name"] == batch_name)
        mask2 = (adata.obs["cell_type"] == cell_types[1]) & (adata.obs["batch_name"] == batch_name)
    else:
        mask1 = (adata.obs["cell_type"] == cell_types[0])
        mask2 = (adata.obs["cell_type"] == cell_types[1])
    return np.array(mask1), np.array(mask2)
'''

def get_groups_masks(adata, cell_types, batch_names = None):
    assert len(cell_types) == 2, "two cell types should be specified."
    if batch_names is not None:
        assert len(batch_names) == 2, "two batch names should be specified."
    if batch_names is not None:
        mask1 = (adata.obs["cell_type"] == cell_types[0]) & (adata.obs["batch_name"] == batch_names[0])
        mask2 = (adata.obs["cell_type"] == cell_types[1]) & (adata.obs["batch_name"] == batch_names[1])
    else:
        mask1 = (adata.obs["cell_type"] == cell_types[0])
        mask2 = (adata.obs["cell_type"] == cell_types[1])
    return np.array(mask1), np.array(mask2)


def normalize(adata, norm_type):
    
    if norm_type == "raw":
        norm_layer_name = "X"
    elif norm_type == "dist":
        norm_counts = np.zeros(adata.shape)
        for gap in (adata.var["bin2_id"] - adata.var["bin1_id"]).unique():
            gap_idx = np.where(adata.var["bin2_id"] - adata.var["bin1_id"] == gap)[0]
            curr_counts = adata.X[:, gap_idx].copy()
            norm_counts[:, gap_idx] = qnorm.quantile_normalize(curr_counts.toarray(), axis = 0)
        norm_layer_name = "norm"
        adata.layers[norm_layer_name] = norm_counts
        
    elif norm_type == "total":
        norm_counts = qnorm.quantile_normalize(adata.X.toarray(), axis = 0)
        norm_layer_name = "norm"
        adata.layers[norm_layer_name] = norm_counts
    else:
        raise ValueError("invalid normalization.")
    return adata, norm_layer_name

def read_anndata_filelist(filelist_path):
    all_anndata_paths = {}
    with open(filelist_path, "r") as f:
        for line in f:
            anndata_filepath, chrom = line.rstrip("\n").split("\t")
            if anndata_filepath in all_anndata_paths.keys():
                all_anndata_paths[anndata_filepath].append(chrom)
            else:
                all_anndata_paths[anndata_filepath] = [chrom]
    return all_anndata_paths

def create_parser():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-f', '--filelist', action = 'store', required = True, 
                        help = 'input filelist which is a two column tab-sep file. First column: anndata path, second column: chromosome name existing in anndata.')
    
    parser.add_argument('-a', '--groupA', action = 'store', required = True,
                        help = 'group A to be compared.')
    
    parser.add_argument('-b', '--groupB', action = 'store', required = True,
                        help = 'group B to be compared.')
    
    parser.add_argument('-s', '--downsample', action = 'store', type = int,
                        help = "number of samples to downsample.")
    
    parser.add_argument('-n', '--norm', action = 'store',
                        help = "normalization type can be raw, total or dist.",
                        default = "total")
    
    parser.add_argument('--batch', action = 'store',
                        help = "if it is set, cell types from a specified batch are compared.")
    
    parser.add_argument('--batch2', action = 'store',
                        help = "if it is set, this batch is used for selection of cells from the second cell type.")
    
    parser.add_argument('-o', '--output-path', action = 'store',
                        help = "output path to store t-tests results.")
    
    return parser


if __name__ == "__main__":
    main()
