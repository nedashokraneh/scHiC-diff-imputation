# this script is based on https://github.com/yezhengSTAT/scVI-3D/tree/master

import os
import json
import shutil
import argparse
from tqdm import tqdm

import scvi
import cooler
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import lil_matrix


def main():

    parser = create_parser()
    args = parser.parse_args()
    config = read_config(args.config_path)
    
    if not os.path.exists(config["output_dir"]):
        os.mkdir(config["output_dir"])
        
    chrom_sizes = read_chrom_lens(config["chrom_size_filepath"])
    adata = ad.read_h5ad(config["adata_path"])
    
    if not all( [ch in adata.var["chrom"].unique() for ch in config["chromosomes"]] ):
        raise ValueError("Not all required chromosomes exist in adata.")
    
    
    
    for chrom in config["chromosomes"]:
        
        idx = np.where(adata.var["chrom"] == chrom)[0]
        adata_ = adata[:, idx].copy()
        imputed_adata_X = lil_matrix(adata_.shape)
        latents = {}
    
        if config["flattening"] == "chrom":
            
            latent, imputed_counts = train(adata_, config)
            imputed_adata_X[:, :] = imputed_counts
            latents["whole"] = latent
        
        elif config["flattening"] == "band":
            
            bands_lists = get_bands_lists(chrom, chrom_sizes, config)
            for band_id, band_num in enumerate(bands_lists):
                band_idx = np.where((adata_.var["bin2_id"] - adata_.var["bin1_id"]) == band_num)[0]
                band_adata = adata_[:, band_idx].copy()
                latent, imputed_counts = train(band_adata, config)
                imputed_adata_X[:, band_idx] = imputed_counts
                latents["band{}".format(band_id + 1)] = latent
            
        else:
            
            bands_lists = get_bands_lists(chrom, chrom_sizes, config)
            for pool_id, band_bounds in enumerate(bands_lists):
                min_band, max_band = band_bounds["min_band"], band_bounds["max_band"]
                pool_idx = np.where(((adata_.var["bin2_id"] - adata_.var["bin1_id"]) >= min_band) &
                                    ((adata_.var["bin2_id"] - adata_.var["bin1_id"]) <= max_band))[0]
                pool_adata = adata_[:, pool_idx].copy()
                print(pool_adata.shape)
                latent, imputed_counts = train(pool_adata, config)
                imputed_adata_X[:, pool_idx] = imputed_counts
                latents["pool{}".format(pool_id + 1)] = latent
                
            
        imputed_adata_X = imputed_adata_X.tocsr()
        imputed_adata = ad.AnnData(imputed_adata_X)
        imputed_adata.obs = adata_.obs
        imputed_adata.var = adata_.var
        for key in latents.keys():
            if latents[key] is not None:
                imputed_adata.obsm[key] = latents[key]
        for k, key in enumerate(imputed_adata.obsm.keys()):
            if k == 0:
                full_latent = imputed_adata.obsm[key]
            else:
                full_latent = np.concatenate([full_latent, 
                                              imputed_adata.obsm[key]], axis = 1)
    
        imputed_adata.obsm["whole"] = full_latent
        imputed_adata_path = os.path.join(config["output_dir"], "{}_scVI.h5ad".format(chrom))
        imputed_adata.write(imputed_adata_path)

def train(adata, config):
    
    orig_adata_shape = adata.shape
    imputed_counts = np.zeros(orig_adata_shape)
    
    valid_cells = list(np.where(adata.X.sum(axis = 1) > 0)[0])
    
    if len(valid_cells) == 0:
    
        latent = None
        
    else:
        
        adata = adata[valid_cells, :].copy()
        valid_feats = list(np.where(adata.X.sum(axis=0) > 0)[1])
        adata = adata[:, valid_feats].copy()
        avg_seq_depth = adata.X.sum(axis = 1).mean()
        
        if config["batch_correction"]:
            scvi.model.SCVI.setup_anndata(adata, batch_key = config["batch_name"])
        else: 
            scvi.model.SCVI.setup_anndata(adata)
            
        model = scvi.model.SCVI(adata, n_latent = config["latent_dim"])
        model.train(use_gpu = config["gpu"])
        
        latent = np.zeros((orig_adata_shape[0], config["latent_dim"]))
        latent[valid_cells, :] = model.get_latent_representation()
        
        if config["batch_correction"]:
            available_batches = list(set(adata.obs[config["batch_name"]]))
            for batch_name in available_batches:
                sub_imputed_counts = model.get_normalized_expression(
                    library_size = avg_seq_depth,
                    transform_batch = batch_name)
                imputed_counts[np.ix_(valid_cells, valid_feats)] += sub_imputed_counts
            imputed_counts /= len(available_batches)
        else:
            sub_imputed_counts = model.get_normalized_expression(
                library_size = avg_seq_depth)
            imputed_counts[np.ix_(valid_cells, valid_feats)] = sub_imputed_counts
    
    return latent, imputed_counts
        

def read_config(config_path):
    
    config = json.load( open(config_path, "r") )
    
    assert config["flattening"] in ["chrom", "band", "pool"], "'chrom', 'band', 'pool'-based flattenings are supported."
    if not "max_distance" in config:
        config["max_distance"] = None
    if not "gpu" in config:
        config["gpu"] = False
    else:
        config["gpu"] = bool(config["gpu"])
    if not "latent_dim" in config:
        config["latent_dim"] = 50
    if not "batch_correction" in config:
        config["batch_correction"] = False
    if config["batch_correction"]:
        if not "batch_name" in config:
            raise ValueError("batch_name should be specified for the batch correction!")
    return config

def read_chrom_lens(f):
    chrom_lens = {}
    with open(f, 'r') as f:
        for line in f:
            chrom, size = line.rstrip('\n').split('\t')
            chrom_lens[chrom] = int(size)
    return chrom_lens

def get_bands_lists(chrom, chromosome_sizes, config):
    
    bands_lists = []
    chrom_num_bands = int(np.ceil(chromosome_sizes[chrom] / config["resolution"]))
    if config["max_distance"] == None:
        max_band = chrom_num_bands
    else:
        max_band = int(config["max_distance"] / config["resolution"])
        max_band = np.minimum(max_band, chrom_num_bands)
        
        
    if config["flattening"] == "band":
        
        for b in np.arange(1, max_band, dtype = "int"):
            bands_lists.append(b)
    
    else:
        
        curr_pool_id = 1
        curr_pool_size = 1
        curr_start_band = 1
        
        while curr_start_band <= max_band:
            #bands_lists.append(np.arange(curr_start_band, curr_start_band + curr_pool_size, dtype = "int"))
            bands_lists.append({'min_band': curr_start_band, 
                                'max_band': curr_start_band + curr_pool_size - 1})
            curr_start_band += curr_pool_size
            curr_pool_id += 1
            curr_pool_size += 1
                
    return bands_lists

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config-path', action = 'store', required = True, \
                        help = 'config path')
    return parser


if __name__ == "__main__":
    main()

