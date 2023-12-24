# this script is based on https://github.com/HuMingLab/SnapHiC/blob/master/src/rwr.py

import os
import gc
import sys
import time
import argparse
from mpi4py import MPI

import re
import cooler
import numpy as np
import pandas as pd
from tqdm import tqdm
import networkx as nx
from scipy.sparse import eye, coo_matrix, csc_matrix


def main():

    parser = create_parser()
    args = parser.parse_args()

    mpi_comm = MPI.COMM_WORLD
    mpi_size = mpi_comm.Get_size()
    mpi_rank = mpi_comm.Get_rank()

    window_size = args.window_size
    step_size = int(window_size/2)
    chromosome_lengths = read_chrom_lens(args.chrom_lens)

    file_paths = []
    if getattr(args, 'filelist') != None:
        with open(args.filelist, 'r') as f:
            for line in f:
                file_paths.append(line.rstrip('\n'))
        file_paths = np.array(file_paths)
    elif getattr(args, 'indir') != None:
        file_paths = [os.path.join(args.indir, fn) for fn in os.listdir(args.indir)]
        file_paths = np.array(file_paths)
    else:
        raise ValueError("One of filelist or indir arguments should be set.")

    if (args.format != 'tab-sep') & (args.format != 'cooler'):
        raise ValueError('cooler and tab-sep input formats are available.')

    if args.format == 'tab-sep':
        if (getattr(args, 'chrom_columns') == None) | (getattr(args, 'pos_columns') == None):
            raise ValueError('chrom-columns and pos-columns are required if input format is tab-sep.')

    tasks = balanced_split(file_paths, mpi_size)
    task_filepaths = tasks[mpi_rank]
    #print(task_filepaths)
    if (mpi_rank == 0) & (not os.path.exists(args.output_dir)):
        os.mkdir(args.output_dir)
    mpi_comm.barrier()


    print('{} tasks assigned to rank {}'.format(len(task_filepaths), mpi_rank))

    for filepath in task_filepaths:
        process_cell(filepath, args.format, args.extension, 
                     args.output_dir, args.resolution, 
                     chromosome_lengths, 
                     args.chrom_columns, args.pos_columns, 
                     args.partial, window_size, step_size, 
                     args.keep_short_range, args.upper_distance, 
                     args.rp)
    mpi_comm.barrier()


def process_cell(filepath, format, extension, output_dir, resolution, chromosome_lengths, chrom_columns, pos_columns, partial, window_size, step_size, keep_short_range, upper_distance, rp):

    if filepath.endswith('::resolutions/{}'.format(resolution)):
        filepath_ = filepath.split('::')[0]
    else:
        filepath_ = filepath
    reg = r"(?P<fname>.*)(%s)" % extension
    filename = re.match(reg, os.path.basename(filepath_)).groupdict()['fname']
    #out_filepath = os.path.join(output_dir, '{}_rwr.cool'.format(filename))
    out_filepath = os.path.join(output_dir, '{}.cool'.format(filename))
    if os.path.exists(out_filepath):
        return
    start_time = time.time()
    bins, offset = make_bins(chromosome_lengths, resolution)
    all_imp_edgelist = pd.DataFrame()
    if format == 'cooler':
        c = cooler.Cooler(filepath)
        coarsening_ratio = int(resolution/c.binsize)
        if coarsening_ratio < 1:
            raise ValueError("resolution cannot be less than cooler bin size.")
        for chrom in chromosome_lengths.keys():
            edgelist = c.matrix(balance = False, as_pixels=True).fetch(chrom)
            edgelist.loc[:, ['bin1_id', 'bin2_id']] -= c.offset(chrom)
            if coarsening_ratio != 1:
                edgelist.loc[:, ['bin1_id', 'bin2_id']] = (edgelist.loc[:, ['bin1_id', 'bin2_id']] / coarsening_ratio).astype(int)
                edgelist = edgelist.groupby(['bin1_id', 'bin2_id']).agg({'count': 'sum'}).reset_index()
            imp_edgelist = rwr(edgelist, chromosome_lengths[chrom], 
                               resolution, 
                               partial, window_size, step_size, 
                               keep_short_range, upper_distance, 
                               rp)
            imp_edgelist.iloc[:, [0,1]] += offset[chrom]
            all_imp_edgelist = pd.concat([all_imp_edgelist, 
                                          imp_edgelist], 
                                         axis = 0)

    elif format == 'tab-sep':
        all_edgelist = pd.read_csv(filepath, sep = "\t", header = None)
        for chrom in chromosome_lengths.keys():
            edgelist = all_edgelist[(all_edgelist.iloc[:, chrom_columns[0]] == chrom) & (all_edgelist.iloc[:, chrom_columns[1]] == chrom)]
            edgelist = edgelist.iloc[:, pos_columns]
            edgelist = (edgelist / resolution).astype(int)
            edgelist.columns = ['bin1_id', 'bin2_id']
            edgelist = edgelist.groupby(['bin1_id', 'bin2_id']).size().to_frame('count').reset_index()
            imp_edgelist = rwr(edgelist, chromosome_lengths[chrom], 
                               resolution, 
                               partial, window_size, step_size, 
                               keep_short_range, upper_distance, 
                               rp)
            imp_edgelist.iloc[:, [0,1]] += offset[chrom]
            all_imp_edgelist = pd.concat([all_imp_edgelist, 
                                          imp_edgelist], 
                                         axis = 0)



    #print('finished imputation of {}.\n'.format(filepath))
    
    cooler.create_cooler(out_filepath,
                    bins=bins,
                    pixels=all_imp_edgelist,
                    ordered=True,
                    dtypes={'count': np.float64})
    print("--- it took %s seconds to process one single cell ---\n" 
          % (time.time() - start_time))
    return


def rwr(edges, chrom_len, resolution, 
        partial, window_size, step_size, 
        keep_short_range, distance_threshold, rp):

    edges.columns = ['bin1_id', 'bin2_id', 'count']
    valid_bins = np.unique(list(edges['bin1_id'].values) + 
                           list(edges['bin2_id'].values))
    chrom_size = int(np.ceil(chrom_len/resolution))
    ws, ss = int(window_size/resolution), int(step_size/resolution)
    distance_threshold = int(distance_threshold/resolution)
    #r = coo_matrix(((0,), ((0,), (0,))), shape = (chrom_size,chrom_size)).todense()
    rwr_edges = pd.DataFrame()
    if partial:
        window_starts = range(0, chrom_size, ss)
    else:
        window_starts = [0]
        #ws = chrom_len
        ws = chrom_size
    for window_start in tqdm(window_starts):
        window_end = window_start + ws
        window_edges = edges[(edges['bin1_id'] >= window_start) & 
                             (edges['bin2_id'] < window_end)]
        
        neighbor_edges = pd.DataFrame({
            'bin1_id': list(range(window_start, 
                                  min(chrom_size-1, window_end - 1))), 
            'bin2_id': list(range(window_start+1, 
                                  min(chrom_size, window_end)))})
        
        window_edges = pd.concat([neighbor_edges[['bin1_id', 'bin2_id']], 
                                  window_edges[['bin1_id', 'bin2_id']]], 
                                 axis = 0)
        
        if window_edges.shape[0] == 0:
            continue
        window_edges.loc[:,'count'] = 1
        
        g = get_stochastic_matrix_from_edgelist(window_edges)
        window_edges_imp = solve_rwr_iterative(g, rp)
        #if partial:
        window_edges_imp = window_edges_imp[(window_edges_imp['i'] + window_edges_imp['j'] > ss) & (window_edges_imp['i'] + window_edges_imp['j'] < ws + ss)]
        if keep_short_range:
            window_edges_imp = window_edges_imp[(window_edges_imp['j'] - window_edges_imp['i'] <= distance_threshold)]
        window_edges_imp['i'] += window_start
        window_edges_imp['j'] += window_start
        window_edges_imp.columns = ['bin1_id', 'bin2_id', 'count']
        window_edges_imp = window_edges_imp[window_edges_imp['bin2_id'] >= window_edges_imp['bin1_id']]
        rwr_edges = pd.concat([rwr_edges, window_edges_imp], axis = 0)
        #partial_r = coo_matrix((window_edges_imp['v'], (window_edges_imp['i'], window_edges_imp['j'])), shape = (chrom_size, chrom_size))
        #r += partial_r
    #r = coo_matrix(r)
    #rwr_edges = pd.DataFrame({'bin1_id':r.row, 'bin2_id': r.col, 'count': r.data})
    #rwr_edges = rwr_edges[rwr_edges['bin2_id'] >= rwr_edges['bin1_id']]
    rwr_edges = rwr_edges[(rwr_edges['bin1_id'].isin(valid_bins)) & (rwr_edges['bin2_id'].isin(valid_bins))]
    return rwr_edges

def get_stochastic_matrix_from_edgelist(edgelist):
    g = nx.from_pandas_edgelist(edgelist, 
                                source = 'bin1_id', 
                                target = 'bin2_id', 
                                edge_attr = ['count'], 
                                create_using = nx.Graph())
    
    degrees = np.array([g.degree(i) for i in g.nodes()])
    m = csc_matrix(nx.adjacency_matrix(g).astype(float))
    m.data /= degrees[m.indices] #stochastic matrix
    del g, degrees
    return m

def solve_rwr_iterative(stoch_matrix, alpha = 0.05, max_iter = None):
    max_iter = max_iter if max_iter else float("inf")
    # what is this for?
    # gc.collect()
    I = eye(stoch_matrix.shape[0], dtype=np.float32)
    delta = float("inf")
    A = I.copy()
    counter = 0
    while delta > 1e-6 and counter < max_iter:
        counter += 1
        Aold = A.copy()
        A = (1-alpha) * stoch_matrix * Aold + (alpha) * I
        delta = (abs(A - Aold)).max()

    A += A.transpose()
    A = coo_matrix(A)
    df = pd.DataFrame({'i':A.row, 'j': A.col, 'count': A.data})
    del A
    return df


def read_chrom_lens(f):
    chrom_lens = {}
    with open(f, 'r') as f:
        for line in f:
            chrom, size = line.rstrip('\n').split('\t')
            chrom_lens[chrom] = int(size)
    return chrom_lens

def make_bins(chromosome_lengths, resolution):
    bins = pd.DataFrame()
    offset = {}
    for chrom in chromosome_lengths.keys():
        offset[chrom] = bins.shape[0]
        starts = np.arange(0,chromosome_lengths[chrom],resolution)
        ends = starts[1:]
        ends = np.insert(ends, ends.shape, chromosome_lengths[chrom])
        bins = pd.concat([bins, pd.DataFrame({'chrom': chrom, 'start': starts, 'end': ends})], axis = 0)
    bins.index = range(bins.shape[0])
    return bins, offset

def balanced_split(lst, n):

    nn = int(len(lst) / n)
    extra = len(lst) - nn*n
    sublsts_lengths = [nn+1] * extra + [nn] * (n-extra)
    sublsts = [lst[start:end] for start, end in zip(np.cumsum([0] + sublsts_lengths[:-1]), np.cumsum(sublsts_lengths))]
    return sublsts

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filelist', action = 'store', required = False, \
                        help = 'path of the file including filepaths to impute')
    parser.add_argument('--indir', action = 'store', required = False, \
                        help = 'path of the folder including raw files')
    parser.add_argument('-o', '--output-dir', action = 'store', required = True, help = 'output directory to save imputed contact maps')
    parser.add_argument('-r', '--resolution', action = 'store', required = True, type = int, help = "resolution of contact maps")
    parser.add_argument('-l', '--chrom-lens', action = 'store', required = True, help = 'path to the chromosome lengths file')
    parser.add_argument('--format', action = 'store', required = True, help = 'format of input files: cooler or tab-sep')
    parser.add_argument('--partial', action = 'store_true', default = False, help = 'partial RWR')
    parser.add_argument('--window-size', action = 'store', default = 10000000, type = int, help = 'window size for window-slicing rwr')
    parser.add_argument('--rp', action = 'store', default = 0.05, type = float, help = 'restart probability')
    parser.add_argument('--chrom-columns', action = 'store', required = False, nargs = 2, type = int, help = 'two integer column numbers for chromosomes')
    parser.add_argument('--pos-columns', action = 'store', required = False, nargs = 2, type = int, help = 'two integer column numbers for positions')
    parser.add_argument('--extension', action = 'store', required = True, help = 'extension of input files like .cool and .txt.gz')
    parser.add_argument('--keep-short-range', action = 'store_true', default = False, \
                        help = 'if set, only saves short range interactions within upper distance', required = False)
    parser.add_argument('--upper-distance', action = 'store', default = 2000000, type = int, help = 'maximum distance between bin pairs to store')

    return parser


if __name__ == "__main__":
    main()
