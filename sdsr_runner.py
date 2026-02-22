from sdsr_inner_functions import *
import os


stdr_instant = STDR_for_Multiple_Alignments(
    sequences_path = os.path.join('alignments', 'all-genes.phylip'), # right now it is 50 genes
    format_type='phylip_multiple_genes', # one phylip file with multiple MSAs in it
    number_of_genes=3, #define the number of gene we want to work with
    min_percent_for_balanced_kmean=0.01, #the parameter that control the balanceness of the partitions.
    inner_method='astral', # 'astral' or 'caml'
    threshold=25, # the parameter that decide when to stop the recursion
    sum_or_median='sum', 
    path_to_exported_trees='sdsr_tree.newick',
    output_estimated_ml_trees_path=os.path.join("estimated_ml_trees", "cur_experiment" + ".newick"), 
    output_astral_file_name=os.path.join("astral_files", "cur_experiment"),
    cd_path_for_astral_exe=os.path.abspath("./spectraltree/libs/astral_bin"), 
    different_num_of_genes_to_reconstruct=True,
    number_of_genes_for_reconstruction=2, # or None
    partition_only = False,
    astral_bin = "astral4"  
)

estimated_tree, total_time_pre_processing, total_time_create_similarity, total_time_partitioning, total_time_merging, total_time_recunstructing = stdr_instant.deep_spectral_tree_reconstruction()
time_profile = f'pre_processing_time: {total_time_pre_processing}, total_time_create_similarity: {total_time_create_similarity}, total_time_partitioning: {total_time_partitioning}, total_time_merging: {total_time_merging}, recunstruction_time: {total_time_recunstructing}'
print(f'estimated_tree: {estimated_tree}')
print(time_profile)