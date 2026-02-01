
from sdsr_inner_function import *
from Bio import AlignIO
from Bio import Phylo
from io import StringIO
import spectraltree
import numpy as np
import dendropy
from k_means_constrained import KMeansConstrained
from more_itertools import locate
import time
import os
from dendropy import Tree, TaxonNamespace
import ete3
import subprocess
import psutil
from concurrent.futures import ProcessPoolExecutor, as_completed
import threading
import datetime
import platform
from sklearn.metrics.cluster import rand_score
import sys
import shutil
import copy  
import time
from sdsr_inner_function import *
from sdsr_inner_function import STDR_for_Multiple_Alignments

start_time = time.time()


# 50 genes
type_of_run = 'pc'
sequences_path = os.path.join('alignments', 'all-genes.phylip')
format_type = 'phylip_multiple_genes' # one phylip file with multiple MSAs in it
number_of_genes = 2 #define the number of gene we want to work with
min_percent_for_balanced_kmean = 0.01 # the parameter that control the balanceness of the partitions.
inner_method = 'caml' # or 'astral'
threshold = 25 # the parameter that decide when to stop the recursion
experiment_name = f'sdsr_with_{inner_method}_from_{number_of_genes}_genes'
path_to_exported_trees = os.path.join('final_tree', experiment_name +'.newick')
output_estimated_ml_trees_path = os.path.join("estimated_ml_trees", experiment_name + ".newick")
output_astral_file_name = os.path.join("astral_files", experiment_name)
output_astral_file_name_linux = output_astral_file_name
output_estimated_ml_trees_path_linux = output_estimated_ml_trees_path
cd_path_for_astral_exe = "/mnt/c/Users/ortalrex/'OneDrive - Intel Corporation'/Desktop/Ortal/Thesis/for_species/New_astral/ASTER-Linux"

stdr_instant = STDR_for_Multiple_Alignments(sequences_path, format_type = format_type, number_of_genes = number_of_genes, 
                                            min_percent_for_balanced_kmean = min_percent_for_balanced_kmean, inner_method = inner_method, 
                                            threshold = threshold, sum_or_median = 'sum', path_to_exported_trees = path_to_exported_trees,
                                            output_estimated_ml_trees_path = output_estimated_ml_trees_path, output_astral_file_name = output_astral_file_name,
                                            cd_path_for_astral_exe = cd_path_for_astral_exe, 
                                            output_estimated_ml_trees_path_linux = output_estimated_ml_trees_path_linux, 
                                            output_astral_file_name_linux = output_astral_file_name_linux,
                                            type_of_run = type_of_run)

estimated_tree, total_time_pre_processing, total_time_create_similarity, total_time_partitioning, total_time_merging, total_time_recunstructing = stdr_instant.deep_spectral_tree_reconstruction()
time_profile = f'pre_processing_time: {total_time_pre_processing}, total_time_create_similarity: {total_time_create_similarity}, total_time_partitioning: {total_time_partitioning}, total_time_merging: {total_time_merging}, recunstruction_time: {total_time_recunstructing}'
