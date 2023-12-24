import os
import gc
import json
import shutil
import argparse
from tqdm import tqdm
from mpi4py import MPI

import cooler
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import lil_matrix, csr_matrix

def main():
    
    parser = create_parser()
    args = parser.parse_args()
    mpi_comm = MPI.COMM_WORLD
    mpi_size = mpi_comm.Get_size()
    mpi_rank = mpi_comm.Get_rank()
    
    config = read_config(args.config)
    filter_regions = read_filter_regions(config["filter_regions_filepath"], 
                                         config["filter_file_resolution"],
                                         config["resolution"])
    sc_filepaths, cell_features = read_summary_file(config["summary_file_path"])
    
    
    coordinates, chrom_valid_bins, chrom_bin_offsets, chrom_lp_offsets, num_lps = make_coordinates(
        config["chrom_size_filepath"], config["chromosomes"], 
        filter_regions, config["resolution"], config["max_distance"])
    
    splitted_cell_ids = balanced_split(np.arange(len(sc_filepaths)), mpi_size)
    my_cell_ids = splitted_cell_ids[mpi_rank]
    
    print("{} cells assigned to rank {}.".format(len(my_cell_ids), mpi_rank))
    if (mpi_rank == 0) & (not os.path.exists(os.path.join(config["output_dir"], 
                                       "temp"))):
        os.mkdir(os.path.join(config["output_dir"], "temp"))
        
    mpi_comm.barrier()
    
    temp_outpath = os.path.join(config["output_dir"], 
                                "temp", "{}.h5ad".format(mpi_rank))
    
    if not os.path.exists(temp_outpath):
        #count_mat = lil_matrix((len(my_cell_ids), num_lps))
        cell_idx = []
        cell_vars_idx = []
        cell_counts = []
        max_gap = int(config["max_distance"] / config["resolution"]) if config["max_distance"] != None else None
        
        
        
        for cell_index, cell_id in tqdm(enumerate(my_cell_ids)):
            
            sc_filepath = sc_filepaths[cell_id]
            if config['format'] == 'tab':
                df = pd.read_csv(sc_filepath, 
                                sep = "\t", header = None, 
                                names = config["colnames"])
                df = df[df["chrom1"] == df["chrom2"]]
                df.drop(columns = ["chrom2"], inplace = True)
                df.rename(columns = {'chrom1': 'chr'}, inplace = True)
                
            elif config['format'] == 'mcool':
                c = cooler.Cooler('{}::/resolutions/{}'.format(sc_filepath, 
                                                                config['resolution']))
            elif config['format'] == 'cool':
                c = cooler.Cooler(sc_filepath)
                
            for chrom in config["chromosomes"]:
                if config['format'] == 'tab':
                    sub_df = df[df["chr"] == chrom].copy()
                    if sub_df.empty:
                        continue
                    sub_df[['bin1_id', 'bin2_id']] = (sub_df[['pos1', 'pos2']] / config["resolution"]).astype(int)
                    first_bin_id, second_bin_id = np.minimum(sub_df['bin1_id'], sub_df['bin2_id']), np.maximum(sub_df['bin1_id'], sub_df['bin2_id'])
                    sub_df['bin1_id'] = first_bin_id
                    sub_df['bin2_id'] = second_bin_id
                    if not 'count' in sub_df.columns:
                        sub_df['count'] = 1
                    sub_df = sub_df.groupby(['chrom1', 'bin1_id', 'bin2_id']).aggregate({'count': 'sum'}).reset_index()
                
                elif config['format'] in ['cool', 'mcool']:
                    try:
                        sub_df = c.matrix(balance = False, as_pixels = True).fetch(chrom)
                    except:
                        continue
                    if sub_df.empty:
                        continue
                    sub_df[['bin1_id', 'bin2_id']] -= c.offset(chrom)
                if max_gap is not None:
                    valid_lps = [(b2 <= (b1 + max_gap)) for (b1, b2) in zip(sub_df['bin1_id'], sub_df['bin2_id'])]
                    sub_df = sub_df[valid_lps]
                valid_lps = [(chrom_valid_bins[chrom][b1] & chrom_valid_bins[chrom][b2]) for (b1, b2) in zip(sub_df['bin1_id'], sub_df['bin2_id'])]
                sub_df = sub_df[valid_lps]
                sub_df['var_idx'] = [(chrom_bin_offsets[chrom][b2] - chrom_bin_offsets[chrom][b1]) + chrom_lp_offsets[chrom][b1] for (b1, b2) in zip(sub_df["bin1_id"], sub_df["bin2_id"])]
                cell_idx.extend([cell_index] * sub_df.shape[0])
                cell_vars_idx.extend(sub_df['var_idx'])
                cell_counts.extend(sub_df['count'])
                #count_mat[cell_index, sub_df['var_idx']] = sub_df['count']
                del sub_df, valid_lps 
            collected = gc.collect()
 
            #print("Garbage collector: collected", "%d objects." % collected)
        
        #count_mat = count_mat.tocsr()
        count_mat = csr_matrix((cell_counts, (cell_idx, cell_vars_idx)), 
                               shape = (len(my_cell_ids), num_lps))
        adata = ad.AnnData(count_mat)
        if len(cell_features) != 0:
            for feat_name in cell_features.keys():
                adata.obs[feat_name] = np.array(cell_features[feat_name])[my_cell_ids]
        adata.obs.index = my_cell_ids
        for var in ['chrom', 'bin1_id', 'bin2_id']:
            adata.var[var] = coordinates[var]
        adata.write(temp_outpath)
        
    
    mpi_comm.barrier()
    
    # TODO: loading all stored anndata with memory for one cpu when n_task > 1
    if mpi_rank == 0:
        
        if (mpi_size == 1) and config["output_format"] == "whole":
            outpath = os.path.join(config["output_dir"], "{}.h5ad".format(config["output_name_prefix"]))
            os.rename(temp_outpath, outpath)
        else:
            
            all_adata = []
            for i in range(mpi_size):
                adata_path = os.path.join(config["output_dir"], "temp", "{}.h5ad".format(i))
                all_adata.append(ad.read_h5ad(adata_path))
            whole_adata = ad.concat(all_adata)
            whole_adata.var = all_adata[0].var
            if config["output_format"] == "whole":
                outpath = os.path.join(config["output_dir"], "{}.h5ad".format(config["output_name_prefix"]))
                whole_adata.write(outpath)
            
            else:
                for chrom in config["chromosomes"]:
                    sub_adata = whole_adata[:, whole_adata.var["chrom"] == chrom].copy()
                    outpath = os.path.join(config["output_dir"], "{}_{}.h5ad".format(config["output_name_prefix"], chrom))
                    sub_adata.write(outpath)

        shutil.rmtree(os.path.join(config["output_dir"], "temp"))

def read_config(config_path):
        
    config = json.load( open(config_path, 'r') )
    if not "max_distance" in config:
        config["max_distance"] = None
    if config["format"] == "tab" and not "colnames" in config:
        raise ValueError("Colnames should be specified for tab format.")
    if not "output_format" in config:
        config["output_format"] = "whole"
    else:
        if not config["output_format"] in ["chrom", "whole"]:
            raise ValueError("output_format should be chrom or whole.")
    return config
    
        
def read_filter_regions(filter_regions_filepath, 
                        filter_file_resolution,
                        resolution):
    
    filter_regions = pd.read_csv(filter_regions_filepath,
                                 sep = "\t", header = None).iloc[:,:2]
    filter_regions.columns = ["chrom", "pos"]
    filter_regions["bin_id"] = (filter_regions["pos"] / resolution).astype(int)
    filter_regions.drop(columns = ["pos"], inplace = True)
    filter_regions = filter_regions.groupby(["chrom", "bin_id"]).size().reset_index()
    filter_regions = filter_regions[filter_regions[0] >= int(resolution / filter_file_resolution)]
    filter_regions.drop(columns = [0], inplace = True)
    #filter_regions.drop_duplicates(inplace = True)
    filter_regions["include"] = False
    return filter_regions
        
    
def make_coordinates(chrom_size_filepath, chromosomes, 
                     filter_regions, 
                     resolution, max_distance):
    
    max_gap = int(max_distance / resolution) if max_distance != None else None
    coordinates = {'chrom': [], 'bin1_id': [], 'bin2_id': []}
    chrom_valid_bins = {}
    curr_bin_offset = 0
    curr_lp_offset = 0
    chrom_bin_offsets = {}
    chrom_lp_offsets = {}
    
    with open(chrom_size_filepath, "r") as f:
        
        for line in f:
            
            chrom, chrom_size = line.rstrip("\n").split("\t")
            chrom_size = int(chrom_size)
            
            if chrom in chromosomes:
                
                starts = np.arange(0, chrom_size, resolution)
                ends = starts + resolution
                ends[-1] = chrom_size
                bin_ids = (starts / resolution).astype(int)
                bins = pd.DataFrame({'chrom': chrom, 
                                     'start': starts, 'end': ends, 
                                     'bin_id': bin_ids})
                bins = pd.merge(bins, filter_regions, on = ["chrom", "bin_id"], how = "left")
                bins.fillna(True, inplace = True)
                num_bins = bins.shape[0]
                chrom_valid_bins[chrom] = np.array(bins["include"])
                
                # bin offsets
                chrom_bin_offsets[chrom] = curr_bin_offset + np.cumsum(bins["include"])
                curr_bin_offset += bins["include"].sum()
                
                # locus-pair offsets
                if max_gap is None:
                    num_valid_lps = [bins["include"][idx] * 
                                     bins["include"][idx : num_bins].sum() 
                                     for idx in range(num_bins)]
                else:
                    num_valid_lps = [bins["include"][idx] * 
                                     bins["include"][idx : (min(idx + max_gap + 1, num_bins))].sum() 
                                     for idx in range(num_bins)]
                #chrom_lp_offsets[chrom] = [0] + list(np.cumsum(num_valid_lps)[ : -1])
                chrom_lp_offsets[chrom] = [curr_lp_offset] + list(np.cumsum(num_valid_lps)[ : -1] + curr_lp_offset)
                curr_lp_offset += np.sum(num_valid_lps)
                
                # coordinates
                bin1_ids = list(np.repeat(np.arange(num_bins), num_valid_lps))
                valid_bin_ids = np.where(bins['include'])[0]
                if max_gap is None:
                    bin2_ids = [bin2_id for bin1_id in valid_bin_ids 
                            for bin2_id in range(bin1_id, num_bins) if bins["include"][bin2_id]]
                else:
                    bin2_ids = [bin2_id for bin1_id in valid_bin_ids 
                                for bin2_id in range(bin1_id, min(bin1_id + max_gap + 1, num_bins)) if bins["include"][bin2_id]]
                coordinates['chrom'] = coordinates['chrom'] + [chrom] * len(bin1_ids)
                coordinates['bin1_id'] = coordinates['bin1_id'] + bin1_ids
                coordinates['bin2_id'] = coordinates['bin2_id'] + bin2_ids
                
    num_lps = curr_lp_offset
    return coordinates, chrom_valid_bins, chrom_bin_offsets, chrom_lp_offsets, num_lps

def read_summary_file(summary_file_path):
    
    sc_filepaths = []
    cell_features = {}
    with open(summary_file_path, 'r') as f:
        first_line = f.readline()
        first_line = first_line.rstrip("\n").split("\t")
        if not first_line[0] == 'filepath':
            raise ValueError('First column of a summary file should include filepaths.')
        if len(first_line) > 1:
            for feat in first_line[1:]:
                cell_features[feat] = []
            feat_names = list(cell_features.keys())   
        for line in f:
            line = line.rstrip("\n").split("\t")
            sc_filepaths.append(line[0])
            if len(line) > 1:
                for i, feat_value in enumerate(line[1:]):
                    cell_features[feat_names[i]].append(feat_value)
            
    return sc_filepaths, cell_features   
    
def balanced_split(lst, n):

    nn = int(len(lst) / n)
    extra = len(lst) - nn*n
    sublsts_lengths = [nn+1] * extra + [nn] * (n-extra)
    sublsts = [lst[start:end] for start, end in zip(np.cumsum([0] + sublsts_lengths[:-1]), np.cumsum(sublsts_lengths))]
    return sublsts


def create_parser():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', action = 'store', required = True, 
                        help = 'path of config including parameters')
    return parser


if __name__ == "__main__":
    main()
