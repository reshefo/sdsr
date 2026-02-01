import time
from Bio import AlignIO
import numpy as np

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

recursion_or_plain = 'recursion'
our_method_name = 'caml'
gene_num = 2
data_set_name = '001'
type_of_datasets = '50_species'
model_name = 'model.50.2000000.0.000001.0'
model_number = '1'
different_num_of_genes_to_recunstruct = False

class STDR_for_Multiple_Alignments():
    
    def __init__(self, path_to_seq, format_type, number_of_genes, min_percent_for_balanced_kmean, inner_method, threshold, sum_or_median, path_to_exported_trees, output_estimated_ml_trees_path, output_astral_file_name, cd_path_for_astral_exe, output_estimated_ml_trees_path_linux, output_astral_file_name_linux, type_of_run):
        self.path_to_seq = path_to_seq
        self.format_type = format_type
        self.number_of_genes = number_of_genes
        self.inner_method = inner_method
        self.min_percent_for_balanced_kmean = min_percent_for_balanced_kmean
        self.threshold = threshold
        self.sum_or_median = sum_or_median
        self.path_to_exported_trees = path_to_exported_trees
        self.output_estimated_ml_trees_path = output_estimated_ml_trees_path
        self.output_astral_file_name = output_astral_file_name
        self.cd_path_for_astral_exe = cd_path_for_astral_exe
        self.output_estimated_ml_trees_path_linux = output_estimated_ml_trees_path_linux
        self.output_astral_file_name_linux = output_astral_file_name_linux
        self.type_of_run = type_of_run
        #self.workers_for_trees = workers_for_trees


    '''
    def show_tree_in_jupyter(self, tree):
    
        from ete3.treeview import TreeStyle
        from PIL import Image
        
        %matplotlib inline
    
        ete_tree = ete3.Tree(str(tree)+';', format=1)
    
        # Configure the tree style (optional)
        ts = TreeStyle()
    
        # Render the tree to an image file
        ete_tree.render("tree.png", w=500, tree_style=ts)
    
        # Display the image in the notebook
        image = Image.open("tree.png")
    
        # Resize the image
        width = 500  # Set the desired width
        height = 1500  # Set the desired height
        resized_image = image.resize((width, height))
    
        # Display the resized image
        return(resized_image)
        '''
    

    def export_tree(self, tree):
        ete_tree_to_export = ete3.Tree(str(tree)+';', format=1) # replace lenth with nan or missing lenth in 1, accure in astral4 (we are interesting in the topology)
        with open(self.path_to_exported_trees, "w") as file:
            file.write(ete_tree_to_export.write(format=1))
            
    #process functions
    '''
    input: row data
    output: clean data: 
        self.alignments
        self.leaves_names
        self.taxa_metadata
    '''
        
    def pre_process(self):
        if self.format_type == 'phylip_multiple_genes':
            self.process_multiple_genes()
        if self.format_type == 'phylip_whole_genome':
            self.process_whole_genome()
        if self.format_type == 'phylip_multiple_genes_from_different_files':
            self.process_multiple_genes_from_different_files()

    
    def clean(self, selected_alignments):
        #create leaves names
        alignment_for_leaves_names = selected_alignments[0]
        alignment_for_leaves_names.sort()
        leaves_names = [str(record.id).split()[0] for record in alignment_for_leaves_names]
        
        if self.format_type == 'phylip_multiple_genes_from_different_files':
            leaves_names = leaves_names      # no outgroup removal needed
        else:
            leaves_names = leaves_names[1:] #remove the outgroup name

        #clean the selected alignments and create clean_alignments, taxa_metadata and leaves_names give one alignment for sorted names (leaves_nams)
        clean_alignments = []
    
        for i in list(range(0,len(selected_alignments))):
            alignment = selected_alignments[i]
            #print('alignment befor sort:')
            #print(alignment)
            alignment.sort()
            #print('alignment after sort:')
            #print(alignment)
            
            if self.format_type == 'phylip_multiple_genes_from_different_files':
                alignment = alignment                     # no outgroup removal needed
            else:
                #remove the outgroup
                alignment = alignment[1:]

            alignment = [str(record.seq) for record in alignment]      #convert the alignment to string
            alignment = list(map(lambda el:[el], alignment))           #convert the alignment to lists


            temp_alignment = []                                        #convert each sequence to list
            for i, element in enumerate(alignment):
                o = list(alignment[i][0])                         
                temp_alignment.append(o)
            alignment = np.array(temp_alignment)                       #convert the alignment to np.array
    
            clean_alignments.append(alignment)
        
        
        taxa_metadata = spectraltree.TaxaMetadata(dendropy.TaxonNamespace(leaves_names),leaves_names)
        
        self.alignments = clean_alignments
        self.leaves_names = leaves_names
        self.taxa_metadata = taxa_metadata

        del temp_alignment
        del alignment
        del clean_alignments
        del leaves_names
        del taxa_metadata
        
        #print('self.alignments:')
        #print(self.alignments)
        #print('self.alignments[0].shape:')
        #print(self.alignments[0].shape)
        #print('self.leaves_names')
        #print(self.leaves_names)

              
    def process_multiple_genes(self):
    
        #read the phylip file and create a list of all the 1000 sequences
        all_alignments = list(AlignIO.parse(self.path_to_seq, 'phylip'))
        
        #define the number of gene we want to work with
        selected_alignments = all_alignments[:self.number_of_genes]

        #clean and create all the needed objects
        self.clean(selected_alignments)

        del all_alignments
        del selected_alignments
        

    def create_gene_windows_from_whole_genom(self, sequence_path, windows_number):
        
        if self.inner_method == 'caml':
            windows = {i: slice(i * 2500, (i + 1) * 2500) for i in range(windows_number)}
            #windows = {i: slice(i * 500, (i + 1) * 500) for i in range(windows_number)}
        if self.inner_method == 'astral':
            windows = {i: slice(i * 500 + i * 2500, (i + 1) * 500 + i * 2500) for i in range(windows_number)}
            
        with open(sequence_path, "r") as alignment_file:
            # Read the entire alignment
            all_alignments = AlignIO.read(alignment_file, "phylip")
    
            alignments = []
            for i in range(windows_number):
                alignments.append(all_alignments[:, windows[i]])

            
            del all_alignments
    
        return alignments

    
    def process_whole_genome(self):
        #print('do process_whole_genome')

        # create gene windows from caster data (2500 nucluid in each window)
        selected_alignments = self.create_gene_windows_from_whole_genom(self.path_to_seq, self.number_of_genes)

        #clean and create all the needed objects
        self.clean(selected_alignments)

        del selected_alignments
        
        #code for proccesing that computes sequences and metadata
        #sequences = np.array([])
        #metadata  = []
        #self.sequences = sequences
        #self.metadata = metadata

    def process_multiple_genes_from_different_files(self):
        # create names of alignments files
        alignments_names_list = []
        for filename in os.listdir(self.path_to_seq):
            name = os.path.join(self.path_to_seq, filename)
            alignments_names_list.append(name)
    
        #define the number of gene we want to work with
        selected_alignments_names = alignments_names_list[:self.number_of_genes]
        print('selected_alignments_names')
        print(selected_alignments_names)
    
        # Read the alignments and append
        selected_alignments = []
        for alignment_file_name in selected_alignments_names:
            alignment = AlignIO.read(alignment_file_name, "phylip")
            selected_alignments.append(alignment)
        #print('process_multiple_genes_from_different_files')
        #print(selected_alignments)
        #clean and create all the needed objects
        self.clean(selected_alignments)
        #print('end cleaning')

        del selected_alignments


    from concurrent.futures import ProcessPoolExecutor, as_completed

    def chunk_list(self, lst, n_chunks):
        # Splits lst into n_chunks as evenly as possible
        avg = len(lst) // n_chunks
        rem = len(lst) % n_chunks
        chunks = []
        start = 0
        for i in range(n_chunks):
            end = start + avg + (1 if i < rem else 0)
            chunks.append(lst[start:end])
            start = end
        return chunks


    
    def process_chunk(self, chunk):
        results = []
        
        # Track start
        start = time.time()
        start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        cpu_core = self.get_cpu_core()
        process = psutil.Process(os.getpid())
    
        for alignment in chunk:
            W = spectraltree.JC_similarity_matrix(alignment)
            W = W.astype(np.float32)
            results.append(W)
    
        # Track end
        end = time.time()
        end_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        

    
        # Duration in seconds + minutes
        duration_sec = round(end - start, 2)
        duration_min = round(duration_sec / 60, 2)
    
        parallel_info = {
            'cpu_core': cpu_core,
            'start_time': start_time,
            'end_time': end_time,
            'duration_sec': duration_sec,
            'duration_min': duration_min,
            'peak_mem_GB': round(peak_gb, 2),
            'alignments_processed': len(chunk)
        }

        print(parallel_info)

        return results, parallel_info


    def create_similarity_matrices(self):
        self.similarities_array = []
        self.parallel_infos = []

        if self.type_of_run == 'pc':
            for alignment in self.alignments:
                W = spectraltree.JC_similarity_matrix(alignment)
                self.similarities_array.append(W)

        elif self.type_of_run == 'cluster':
            peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            mem_after_append = peak_kb / 1024
            print(f"Memory after append:  {mem_after_append:.2f} MB")
            
            self.max_workers = int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 1))
            print(f"Using {self.max_workers} CPU workers")
            
            alignments = self.alignments
            n_workers = 32 #self.max_workers
            chunks = self.chunk_list(alignments, n_workers)
            
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                futures = [executor.submit(self.process_chunk, chunk) for chunk in chunks]
                for future in as_completed(futures):
                    results, parallel_info = future.result()
                    self.similarities_array.extend(results)
                    self.parallel_infos.append(parallel_info)
            
            peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            mem_after_append = peak_kb / 1024
            print(f"Memory after append:  {mem_after_append:.2f} MB")
            
            self.peak_memory_similarity = mem_after_append
            print('self.peak_memory_similarity:')
            print(self.peak_memory_similarity)

            self.print_summary() 

    
    def get_cpu_core(self):
        try:
            with open("/proc/self/stat", "r") as f:
                fields = f.read().split()
                return int(fields[38])
        except Exception:
            return "Unavailable"
    
    def process_alignment(self, i):
        start = time.time()
        cpu_core = self.get_cpu_core()
        process = psutil.Process(os.getpid())

        if type_of_datasets == '10000_species':
            W = spectraltree.JC_similarity_matrix(self.alignments[i])
            W = W.astype(np.float32)
        else:
            W = spectraltree.JC_similarity_matrix(self.alignments[i])


        parallel_info = {
            'index': i,
            'cpu_core': cpu_core,
            'duration': round(time.time() - start, 2),
            'peak_mem GB': round(peak_gb, 2)
        }
        print(parallel_info)
        print('self.alignments[i].shape', self.alignments[i].shape)
        print('len(self.similarities_array):', len(self.similarities_array))
        
        return W, parallel_info
    
    def print_summary(self):
        print('Similarity matrices computed:', len(self.similarities_array))
        '''print('Parallel info:')
        for info in self.parallel_infos:
            print(info)'''


    
    def find_closest_index(self, target_value, values, indexes):
        '''
        function that used to find the index of the "average" outgroup.
        '''
        closest_index = None
        closest_distance = float('inf')  # Initialize with a large value
    
        for index in indexes:
            current_value = values[index]
            distance = abs(target_value - current_value)
            
            if distance < closest_distance:
                closest_distance = distance
                closest_index = index
    
        return closest_index

    def create_fiedler_from_partial_data(self, node):
        laplacians_array = []        
        #print('self.similarities_array')
        #print(self.similarities_array)
        for similarity in self.similarities_array:
            #select the sub group
            if recursion_or_plain == 'find_partition':
                similarity_selected_species = similarity[node.bitmap][:, node.bitmap]
            else:
                similarity_selected_species = similarity[node.bitmap_with_outgroup][:, node.bitmap_with_outgroup]
            #print('type similarity_selected_species')
            #print(type(similarity_selected_species))
            #print('similarity_selected_species')
            #print(similarity_selected_species)
            
            '''
            # Create a degree matrix
            D_normalized = np.diag(np.sum(np.array(similarity_selected_species), axis=1)**(-(1/2)))
            # Create normalized Laplacian matrix
            L = D_normalized @ similarity_selected_species @ D_normalized
            '''
            # new partition
            d_inv = np.sum(similarity_selected_species,axis = 1)**(-0.5)
            L = np.outer(d_inv, d_inv) * similarity_selected_species            

            
            laplacians_array.append(L)        
        # Create Laplacian matrices
        if self.sum_or_median == 'sum':
            weighted_laplacians = sum(laplacians_array)
        if self.sum_or_median == 'median':
            weighted_laplacians = np.median(np.array(laplacians_array), axis=0)
        #print('weighted_laplacians:')
        #print(weighted_laplacians)
        # find the eigenvalus (e) and eigenvectors (v)
        e, v = np.linalg.eigh(weighted_laplacians)

        del similarity_selected_species
        del d_inv #del D_normalized
        del L
        del weighted_laplacians
        del laplacians_array

        partial_fiedler = v[:,len(v)-2]
        self.partial_fiedler = partial_fiedler
    
    def create_partition_by_balanced_kmean(self):

        #apply balanced k-mean
        fiedler = self.partial_fiedler.reshape(-1, 1)           # Reshape the fiedler array to be 2D       
        # Define the constrained k-means clustering
        clf = KMeansConstrained(
            n_clusters=2,
            size_min = int(max(3, len(fiedler) * self.min_percent_for_balanced_kmean)),
            size_max=None,
            random_state=0)
        clf.fit_predict(fiedler)        
        # Get the 2 clusters labels
        partial_partition = clf.labels_.tolist()

        self.partial_partition = partial_partition
        
        
    def splitTaxa(self, node):
        
        self.create_fiedler_from_partial_data(node)    
        self.create_partition_by_balanced_kmean()

        
        #print('self.partial_partition:')
        #print(self.partial_partition)
        
        # find the leaves name in the sub group
        selected_leaves_names = [num for num, flag in zip(self.leaves_names, node.bitmap_with_outgroup) if flag]
        
        #print('self.leaves_names:')
        #print(self.leaves_names)
        #print('selected_leaves_names:')
        #print(selected_leaves_names)
        
        #create bitmaps without outgroup
        leaves_names_without_outgroup_1 = [name for name, flag in zip(selected_leaves_names, self.partial_partition) if flag == 1]
        leaves_names_without_outgroup_2 = [name for name, flag in zip(selected_leaves_names, self.partial_partition) if flag == 0]
        #print("leaves_names_without_outgroup_1:")
        #print(leaves_names_without_outgroup_1)
        #print("leaves_names_without_outgroup_2:")
        #print(leaves_names_without_outgroup_2)
        #find this leaves name in the complete list of names and create 2 bitmaps
        bitmap1 = [leave_name in leaves_names_without_outgroup_1 for leave_name in self.leaves_names]
        bitmap2 = [leave_name in leaves_names_without_outgroup_2 for leave_name in self.leaves_names]

        #print('bitmap1:')
        #print(bitmap1)
        #print('bitmap2:')
        #print(bitmap2)
        
        #buildind sub array 1 with outgrop
        indexes_1 = list(locate(self.partial_partition, lambda x: x == 1))
        indexes_2 = list(locate(self.partial_partition, lambda x: x == 0))        
        #find the index of the average value in group 1
        fiedler_values_in_group_1 = self.partial_fiedler[indexes_1]
        average_fiedler_values_in_group_1 = np.mean(fiedler_values_in_group_1)
        index_most_similar_value_to_average_in_group_1 = self.find_closest_index(average_fiedler_values_in_group_1, self.partial_fiedler, indexes_1)
        #find the index of the average value in group 2
        fiedler_values_in_group_2 = self.partial_fiedler[indexes_2]
        average_fiedler_values_in_group_2 = np.mean(fiedler_values_in_group_2)
        index_most_similar_value_to_average_in_group_2 = self.find_closest_index(average_fiedler_values_in_group_2, self.partial_fiedler, indexes_2)        
        #insert the average index from the other group
        indexes_1.insert(0, index_most_similar_value_to_average_in_group_2)
        indexes_2.insert(0, index_most_similar_value_to_average_in_group_1)

        #creat bitmaps with outgroup
        #find the leaves name of the indecies1 and 2 (in the sub group)
        leaves_name_from_indecies_1 = [selected_leaves_names[i] for i in indexes_1]
        #print('leaves_name_from_indecies_1 with outgroup')
        #print(leaves_name_from_indecies_1)
        leaves_name_from_indecies_2 = [selected_leaves_names[i] for i in indexes_2]
        #print('leaves_name_from_indecies_2 with outgroup')
        #print(leaves_name_from_indecies_2)

        #define outgroup names
        outgroup_name_1 = leaves_name_from_indecies_1[0]
        outgroup_name_2 = leaves_name_from_indecies_2[0]
        #print('outgroup_name_1')
        #print(outgroup_name_1)
        #print('outgroup_name_2')
        #print(outgroup_name_2)
        
        #find this leaves name in the complete list of names and create 2 bitmaps
        bitmap_with_outgroup_1 = [leave_name in leaves_name_from_indecies_1 for leave_name in self.leaves_names]
        bitmap_with_outgroup_2 = [leave_name in leaves_name_from_indecies_2 for leave_name in self.leaves_names]
        #print('leave_name:')
        #print(self.leaves_names)
        
        
        #print('bitmap_with_outgroup_1:')
        #print(bitmap_with_outgroup_1)
        #print('bitmap_with_outgroup_2:')
        #print(bitmap_with_outgroup_2)
        #print('---------')

        return bitmap1, bitmap2, bitmap_with_outgroup_1, bitmap_with_outgroup_2, outgroup_name_1, outgroup_name_2


    

    
    def splitTaxaForPartitionCheck(self, node):
        
        self.create_fiedler_from_partial_data(node)    
        self.create_partition_by_balanced_kmean()

        
        #print('self.partial_partition:')
        #print(self.partial_partition)
        
        # find the leaves name in the sub group
        selected_leaves_names = [num for num, flag in zip(self.leaves_names, node.bitmap) if flag]
        
        #print('self.leaves_names:')
        #print(self.leaves_names)
        #print('selected_leaves_names:')
        #print(selected_leaves_names)
        
        #create bitmaps without outgroup
        leaves_names_without_outgroup_1 = [name for name, flag in zip(selected_leaves_names, self.partial_partition) if flag == 1]
        leaves_names_without_outgroup_2 = [name for name, flag in zip(selected_leaves_names, self.partial_partition) if flag == 0]
        #print("leaves_names_without_outgroup_1:")
        #print(leaves_names_without_outgroup_1)
        #print("leaves_names_without_outgroup_2:")
        #print(leaves_names_without_outgroup_2)
        #find this leaves name in the complete list of names and create 2 bitmaps
        bitmap1 = [leave_name in leaves_names_without_outgroup_1 for leave_name in self.leaves_names]
        bitmap2 = [leave_name in leaves_names_without_outgroup_2 for leave_name in self.leaves_names]

        #print('bitmap1:')
        #print(bitmap1)
        #print('bitmap2:')
        #print(bitmap2)

        
        #print('bitmap_with_outgroup_1:')
        #print(bitmap_with_outgroup_1)
        #print('bitmap_with_outgroup_2:')
        #print(bitmap_with_outgroup_2)
        #print('---------')

        return bitmap1, bitmap2

    
    #reconstruct methods

    def root_with_outgroup_and_remove_it(self, node):
        mrca = self.sub_tree.mrca(taxon_labels=[node.outgroup_name])
        self.sub_tree.reroot_at_edge(mrca.edge, update_bipartitions=True)        
        
        #print('subtree after reroot:')
        #print(self.sub_tree.as_ascii_plot())
        
        self.sub_tree.prune_taxa_with_labels([node.outgroup_name])
        #print('subtree after reroot and remove outgroup:')
        #print(self.sub_tree.as_ascii_plot())
            
    def apply_caml(self, node):
    
        # concate
        concated_alignments = np.empty((self.alignments_selected_rows[0].shape[0], 0))
        for alignment in self.alignments_selected_rows:
            concated_alignments = np.append(concated_alignments, alignment, axis = 1)
        print("Shape of concatenated alignments:", concated_alignments.shape)
            

        #print('concated_alignments(sub array):')
        #print(concated_alignments)
        #print('concated_alignments.shape:')
        #print(concated_alignments.shape)
        
        raxml_instant = spectraltree.RAxML()
        sub_tree = raxml_instant(concated_alignments, self.selected_taxa_metadata)
        
        self.sub_tree = sub_tree
        #print('subtree before reroot:')
        #print(self.sub_tree.as_ascii_plot())
        

    def apply_astral_without_parrelal_genes_reconstruction(self, node):
        ml_trees = []
        for alignment in self.alignments_selected_rows:    
            raxml_instant = spectraltree.RAxML()
            sub_ml_tree = raxml_instant(alignment, self.selected_taxa_metadata)
            ml_trees.append(sub_ml_tree)
        #export each tree to the file
        path_ml_trees = self.output_estimated_ml_trees_path + node.tree_id
        with open(path_ml_trees, 'w') as output_file:
            for tree in ml_trees:
                output_file.write(tree.as_string('newick'))

        # run astral4 in linux
        path_ml_trees_linux = self.output_estimated_ml_trees_path_linux + node.tree_id
        output_astral_file_name = self.output_astral_file_name_linux + node.tree_id
        command = fr'cd {self.cd_path_for_astral_exe} && {astral_name} -u 0 -i {path_ml_trees_linux} -o {output_astral_file_name}' #topology only
        #print(command)
        
        if self.type_of_run == 'pc':
            subprocess.run(['bash', '-c', command], shell=True)
        if self.type_of_run == 'cluster':
            subprocess.run(command, shell=True)

        #import the sub species tree
        sub_tree = Tree.get(file=open(self.output_astral_file_name + node.tree_id, 'r'), schema="newick")   
        self.sub_tree = sub_tree
        #return sub_tree


    
    def reconstruct_single_ml_tree(self, alignment, selected_taxa_metadata, output_path, alignment_index=None):
        """Helper function to reconstruct a single ML tree with RAxML and write to file"""
        # Start timing and memory tracking
        start = time.time()
        start_time = datetime.datetime.now()
        process = psutil.Process()
        
        # Get system info
        cpu_core = self.get_cpu_core()
        node_name = os.environ.get("SLURMD_NODENAME", "local")
        
        # Print initial info
        print(f'alignment_index: {alignment_index}', flush=True)
        print(f'cpu_core: {cpu_core}', flush=True)
        print(f'node_name: {node_name}', flush=True)
        print(f'start_time: {start_time}', flush=True)
        
        # Run RAxML
        raxml_instant = spectraltree.RAxML()
        sub_ml_tree = raxml_instant(alignment, self.selected_taxa_metadata)
        print('self.selected_taxa_metadata')
        print(self.selected_taxa_metadata)
        print('alignment data')
        print(alignment) 
        
        # Write this tree to the output file immediately (append mode)
        with open(output_path, 'a') as output_file:
            output_file.write(sub_ml_tree.as_string('newick'))
        
        # End timing and memory tracking
        end_time = datetime.datetime.now()
        peak_mem = process.memory_info().rss
        peak_gb = peak_mem / (1024 ** 3)
        
        # Print final info
        print(f'alignment_index: {alignment_index} - end_time: {end_time}', flush=True)
        print(f'alignment_index: {alignment_index} - duration: {round(time.time() - start, 2)} seconds', flush=True)
        print(f'alignment_index: {alignment_index} - peak_mem GB: {round(peak_gb, 2)}', flush=True)
        
        # Create info dictionary
        info = {
            'alignment_index': alignment_index,
            'cpu_core': str(cpu_core),
            'node': node_name,
            'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
            'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
            'duration': round(time.time() - start, 2),
            'peak_mem GB': round(peak_gb, 2)
        }
        
        print(f'alignment_index: {alignment_index} - Info: {info}', flush=True)
        
        return sub_ml_tree
    
    def apply_astral(self, node):
        # Print overall job info
        tree_id = node.tree_id
        print('node.tree_id')
        print(node.tree_id)
        
        num_alignments = len(self.alignments_selected_rows)
        
        print(f'=== Starting ASTRAL job ===', flush=True)
        print(f'tree_id: {tree_id}', flush=True)
        print(f'num_alignments: {num_alignments}', flush=True)
        print(f'type_of_run: {self.type_of_run}', flush=True)
        
        # Prepare output path
        path_ml_trees = self.output_estimated_ml_trees_path + node.tree_id
        
        # Clear the file if it exists (start fresh)
        open(path_ml_trees, 'w').close()
        
        ml_trees = []
        
        if self.type_of_run == 'cluster':
            # Parallelize RAxML reconstruction and writing on cluster
            cpus_on_node = int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 1))
            print(f'cpus_on_node: {cpus_on_node}', flush=True)
            
            with ProcessPoolExecutor(max_workers=cpus_on_node) as executor:
                futures = []
                for idx, alignment in enumerate(self.alignments_selected_rows):
                    future = executor.submit(self.reconstruct_single_ml_tree, 
                                            alignment, 
                                            self.selected_taxa_metadata,
                                            path_ml_trees,
                                            idx)
                    futures.append(future)
                
                # Collect results
                for future in futures:
                    ml_trees.append(future.result())
        
        elif self.type_of_run == 'pc':
            # Sequential processing on PC
            for idx, alignment in enumerate(self.alignments_selected_rows):
                ml_tree = self.reconstruct_single_ml_tree(alignment, 
                                                     self.selected_taxa_metadata, 
                                                     path_ml_trees,
                                                     idx)
                ml_trees.append(ml_tree)
        
        print(f'=== All ML trees completed, running ASTRAL ===', flush=True)
        
        # Run astral4 in linux (after all ML trees are generated)
        path_ml_trees_linux = self.output_estimated_ml_trees_path_linux + node.tree_id
        output_astral_file_name = self.output_astral_file_name_linux + node.tree_id
        command = fr'cd {self.cd_path_for_astral_exe} && {astral_name} -u 0 -i {path_ml_trees_linux} -o {output_astral_file_name}'
        
        if self.type_of_run == 'cluster':
            subprocess.run(command, shell=True)
        elif self.type_of_run == 'pc':
            subprocess.run(['bash', '-c', command], shell=True)
        
        # Import the sub species tree
        sub_tree = Tree.get(file=open(self.output_astral_file_name + node.tree_id, 'r'), schema="newick")   
        self.sub_tree = sub_tree
        
        print(f'=== ASTRAL completed for tree_id: {tree_id} ===', flush=True)
    
    def apply_selected_reconstruction_method(self, node):
        if self.inner_method == 'caml':
            self.apply_caml(node)
        if self.inner_method == 'astral':
            self.apply_astral(node)
    
    def reconstruct_alg_wrapper(self, node):
        selected_leaves_names = [leave_name for leave_name, flag in zip(self.leaves_names, node.bitmap_with_outgroup) if flag]
        
        selected_taxa_metadata = spectraltree.TaxaMetadata(dendropy.TaxonNamespace(selected_leaves_names),selected_leaves_names)
        alignments_selected_rows = []
        #print('self.alignments:')
        #print(self.alignments)
        import sys
        #print('size:')
        #print(sys.getsizeof(self.alignments))
        total_bytes = sum(arr.nbytes for arr in self.alignments)
        #print(f"Total bytes: {total_bytes}")
        #print(f"Total MB: {total_bytes / 1024**2:.2f} MB")
        
        if different_num_of_genes_to_recunstruct:
            self.alignments = self.alignments[ :number_of_genes_for_reconstruction]
            print(f'partition run on {self.number_of_genes} and recunstruction run on {int(number_of_genes_for_reconstruction)} genes')
        
        for alignment in self.alignments:
            alignment_selected_rows = alignment[node.bitmap_with_outgroup,:]
            alignments_selected_rows.append(alignment_selected_rows)

        self.alignments_selected_rows = alignments_selected_rows
        self.selected_taxa_metadata = selected_taxa_metadata
        self.selected_leaves_names = selected_leaves_names
        #print('self.selected_leaves_names:')
        #print(self.selected_leaves_names)
        #print('selected_taxa_metadata:')
        #print(selected_taxa_metadata)
        #print('alignments_selected_rows:')
        #print(alignments_selected_rows)
        
        self.apply_selected_reconstruction_method(node)
        if self.make_partition == 'no':
            self.sub_tree = self.sub_tree
        else:
            self.root_with_outgroup_and_remove_it(node)
        #print(self.sub_tree.as_ascii_plot())

        
        return self.sub_tree

    
    def mergeTreesLeftRight(self, node):
        united_tree = dendropy.Tree()
        united_tree.seed_node = dendropy.Node()
        
        united_tree.seed_node.add_child(node.left.tree)        
        united_tree.seed_node.add_child(node.right.tree)

        # convert to dendropy tree
        newick_str = str(united_tree)+';'
        united_tree = dendropy.Tree.get(data=newick_str, schema="newick")
        
        return united_tree

    
    def duplicate_taxa(self, tree1, tree2, taxon_label):
        """
        Checks if a specific taxon label appears in both trees.
        """
        # Collect taxon labels from the first tree
        tree1_taxa = {node.taxon.label for node in tree1.nodes() if node.taxon}
        
        # Collect taxon labels from the second tree
        tree2_taxa = {node.taxon.label for node in tree2.nodes() if node.taxon}
        
        # Check if the taxon label appears in both sets
        if taxon_label in tree1_taxa and taxon_label in tree2_taxa:
            return True
        return False

    
    def mergeTreesLeftRight_for_parallel(self, tree_1, tree_2):
        united_tree = dendropy.Tree()
        united_tree.seed_node = dendropy.Node()
        
        #united_tree.seed_node.add_child(tree_1)        
        #united_tree.seed_node.add_child(tree_2)
        united_tree.seed_node.add_child(tree_1.seed_node)
        united_tree.seed_node.add_child(tree_2.seed_node)

        # convert to dendropy tree
        newick_str = str(united_tree)+';'
        united_tree = dendropy.Tree.get(data=newick_str, schema="newick")
    
        return united_tree
    
     
    def reroot_and_prun_for_parallel(self, merged_tree, outgroup_name_of_upper_level):
        # if there is outgropup (and we are not in the final 2 trees) then reroot acording the outgrop and remove outgrop
        if outgroup_name_of_upper_level is not None:
            mrca = merged_tree.mrca(taxon_labels=[outgroup_name_of_upper_level])
            merged_tree.reroot_at_edge(mrca.edge, update_bipartitions=True)        
            merged_tree.prune_taxa_with_labels([outgroup_name_of_upper_level])
        return merged_tree

    
    def merge_trees_for_parallel(self, nodes_data):
    
        while True:
            # Find all keys longer than 1 character (non-root nodes)
            complex_keys = [key for key in nodes_data if len(key) > 1]
            if not complex_keys:
                break
    
            # Sort keys by length descending (process deepest nodes first)
            complex_keys.sort(key=lambda x: -len(x))
            #print(complex_keys)
            merged = False
    
            for key in complex_keys:
                if key[-1] == 'r':
                    sibling = key[:-1] + 'l'
                elif key[-1] == 'l':
                    sibling = key[:-1] + 'r'
                else:
                    continue
    
                parent = key[:-1]
                #print(parent)
                
                if sibling in nodes_data:
                    # Merge and assign to parent
                    #merged_bitmap = [a or b for a, b in zip(nodes_data[key], nodes_data[sibling])]
                    merged_tree = self.mergeTreesLeftRight_for_parallel(nodes_data[key]['tree'], nodes_data[sibling]['tree'])
                    merged_tree = self.reroot_and_prun_for_parallel(merged_tree, nodes_data[parent]['outgroup_name'])
                    nodes_data[parent]['tree'] = merged_tree
                    #print(nodes_data[parent]['tree'].as_ascii_plot())
                    # Remove merged children
                    del nodes_data[key]
                    del nodes_data[sibling]
    
                    merged = True
                    break  # Restart traversal
    
            if not merged:
                break
    
        # Always set root node if at least one side exists
        if 'r' in nodes_data and 'l' in nodes_data:
            nodes_data['']['tree'] = self.mergeTreesLeftRight_for_parallel(nodes_data['r']['tree'], nodes_data['l']['tree'])
            del nodes_data['r']
            del nodes_data['l']
        elif 'r' in nodes_data:
            nodes_data['']['tree'] = nodes_data['r']['tree']
            del nodes_data['r']
        elif 'l' in nodes_data:
            nodes_data['']['tree'] = nodes_data['l']['tree']
            del nodes_data['l']
        return nodes_data['']['tree']
        
    def create_folder_to_subtrees(self):
        #open a folder to contain the trees
        folder_name_subtrees = (f"{our_method_name}_gen_num_{gene_num}ds_{data_set_name}")
        if type_of_datasets == '50_species':
            folder_name_subtrees = (f"{our_method_name}_{type_of_datasets}_m_{model_number}_{self.min_percent_for_balanced_kmean}_threshold_{self.threshold}_gen_num_{gene_num}ds_{data_set_name}")

        folder_export_subtrees = os.path.join("exported_trees", folder_name_subtrees)
        os.makedirs(folder_export_subtrees, exist_ok=True)
        return folder_export_subtrees
        
    def reconstruct_small_subtrees(self):
        """
        Reconstruct trees for all entries in self.nodes_data where the number of taxa 
        (i.e., True values in 'bitmap') is less than or equal to self.threshold. - ****it can be on single CPU, it written for the structure change
        """
        folder_export_subtrees = self.create_folder_to_subtrees()


        #print('self.nodes_to_reconstruct befor:')
        #print(self.nodes_data)
        for tree_id, data in self.nodes_data.items():
            if np.sum(data["bitmap"]) <= self.threshold:
                node = MyNode(
                    data=data["bitmap"],
                    data_with_outgroup=data["bitmap_with_outgroup"],
                    outgroup_name=data["outgroup_name"],
                    tree_id=tree_id
                )
                reconstructed_tree = self.reconstruct_alg_wrapper(node)

                #save the tree
                tree_file = os.path.join(folder_export_subtrees, f"{tree_id}.nwk")
                print('tree_file path')
                print(tree_file)
                with open(tree_file, 'w') as f:
                    f.write(reconstructed_tree.as_string('newick'))

                #self.nodes_data[tree_id]["tree"] = reconstructed_tree
        
                # print(reconstructed_tree.as_ascii_plot())
        #print('self.nodes_to_reconstruct after:')
        #print(self.nodes_data)


        # retrive back all the trees:
        for tree_id, data in self.nodes_data.items():
            if np.sum(data["bitmap"]) <= self.threshold:
                subtree_file_for_parrelal = os.path.join(folder_export_subtrees, f"{tree_id}.nwk")
                while not os.path.exists(subtree_file_for_parrelal):
                    time.sleep(0.5)  # Sleep half a second before checking again
                # File now exists; load it
                with open(subtree_file_for_parrelal, 'r') as f:
                    loaded_subtree = Tree.get(file=f, schema="newick")
                self.nodes_data[tree_id]["tree"] = loaded_subtree
                
        #print('self.nodes_to_reconstruct after:')
        #print(self.nodes_data)


    
    def _reconstruct_single_tree(self, tree_id, data):

            folder_export_subtrees = self.create_folder_to_subtrees()

            start_time = datetime.datetime.now()
            start = time.time()
            cpu_core = self.get_cpu_core()
            node_name = platform.node()
            bitmap_sum = int(np.sum(data["bitmap"]))

            print(f'tree_id: {tree_id}', flush=True)
            print(f'bitmap_sum: {bitmap_sum}', flush=True)
            print(f'cpu_core: {cpu_core}', flush=True)
            print(f'node_name: {node_name}', flush=True)
            print(f'start_time: {start_time}', flush=True)

            node = MyNode(
                data=data["bitmap"],
                data_with_outgroup=data["bitmap_with_outgroup"],
                outgroup_name=data["outgroup_name"],
                tree_id=tree_id
            )
            reconstructed_tree = self.reconstruct_alg_wrapper(node)
            
            #save the tree
            tree_file = os.path.join(folder_export_subtrees, f"{tree_id}.nwk")
            with open(tree_file, 'w') as f:
                f.write(reconstructed_tree.as_string('newick'))
            print(f'tree of {folder_export_subtrees} id {tree_id} saved')
                    
            end_time = datetime.datetime.now()
            peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            peak_mb = peak_kb / 1024
            peak_gb = peak_mb / 1024

            info = {
                'tree_id': tree_id,
                'bitmap_sum': bitmap_sum,
                'cpu_core': cpu_core,
                'node': node_name,
                'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
                'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                'duration': round(time.time() - start, 2),
                'peak_mem GB': round(peak_gb, 2)
            }
            print('info at the end recunstruction:', flush=True)
            print(info, flush=True)
            print(f"[{tree_id}] About to return result of size: {sys.getsizeof(reconstructed_tree)}", flush=True)
            print(f"at time: {datetime.datetime.now()}", flush=True)

            #return tree_id, reconstructed_tree, info

    def run_astral_for_subtree(self, tree_id, data, folder_export_subtrees):
        """
        Run ASTRAL for a single subtree: read gene trees file, execute ASTRAL, root and remove outgroup.
        Returns the tree_id and the processed tree.
        """
        start = time.time()
        start_time = datetime.datetime.now()
        cpu_core = self.get_cpu_core()
        node_name = platform.node()
        process = psutil.Process()
        
        #print(f'Running ASTRAL for tree_id: {tree_id}', flush=True)
        #print(f'cpu_core: {cpu_core}, node: {node_name}', flush=True)
        
        # Create node for this subtree
        node = MyNode(
            data=data["bitmap"],
            data_with_outgroup=data["bitmap_with_outgroup"],
            outgroup_name=data["outgroup_name"],
            tree_id=tree_id
        )
        
        # Read the gene trees file (already written in Phase 1)
        ml_trees_file = os.path.join(folder_export_subtrees, f"{tree_id}_ml_trees.newick")
        
        # Run ASTRAL
        output_astral_file = os.path.join(folder_export_subtrees, f"{tree_id}_astral_raw.newick")
        
        command = f'cd {self.cd_path_for_astral_exe} && {astral_name} -u 0 -i {ml_trees_file} -o {output_astral_file}'
        
        subprocess.run(command, shell=True)
        
        # Load the ASTRAL tree
        sub_tree = Tree.get(file=open(output_astral_file, 'r'), schema="newick")
        
        # Root with outgroup and remove it
        mrca = sub_tree.mrca(taxon_labels=[node.outgroup_name])
        sub_tree.reroot_at_edge(mrca.edge, update_bipartitions=True)
        sub_tree.prune_taxa_with_labels([node.outgroup_name])
        
        # Track end
        end_time = datetime.datetime.now()
        peak_mem = process.memory_info().rss
        peak_gb = peak_mem / (1024 ** 3)
        
        duration_sec = round(time.time() - start, 2)
        duration_min = round(duration_sec / 60, 2)

        info = {
        'tree_id': tree_id,
        'cpu_core': cpu_core,
        'node': node_name,
        'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
        'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
        'duration_min': duration_min,
        'peak_mem_GB': round(peak_gb, 2)}

        print(info, flush=True)

        return tree_id, sub_tree


    def reconstruct_all_trees_for_single_marker_for_astral(self, marker_idx, trees_in_marker):
        """
        Reconstruct all gene trees for a single marker/alignment.
        Processes all subtrees for this marker sequentially on one CPU.
        
        Args:
            marker_idx: Index of the marker/alignment
            trees_in_marker: List of tasks [(tree_id, data, alignment_idx, alignment), ...]
        
        Returns:
            dict: {tree_id: tree_string, ...} for all subtrees in this marker
        """
        marker_start = time.time()
        cpu_core = self.get_cpu_core()
        node_name = platform.node()
        start_time = datetime.datetime.now()
        
        '''print(f"  Worker started - marker_idx: {marker_idx}, subtrees: {len(trees_in_marker)}, "
              f"cpu_core: {cpu_core}, node: {node_name}, start_time: {start_time}", flush=True)'''
        
        # Dictionary to store results for this marker
        results = {}
        
        # Process each subtree for this marker
        for task_idx, task in enumerate(trees_in_marker):
            tree_id, data, alignment_idx, alignment_selected = task
            
            # Reconstruct the gene tree for this subtree
            _, _, tree_string = self.reconstruct_single_gene_tree_for_subtree_for_astral(
                tree_id,
                data,
                alignment_idx,
                alignment_selected
            )
            
            # Store result
            results[tree_id] = tree_string
                
        
        # Marker processing completed
        marker_end = time.time()
        marker_duration = round((marker_end - marker_start) / 60, 2)
        end_time = datetime.datetime.now()
        
        #successful_trees = sum(1 for v in results.values() if v is not None)
        
        '''print(f"  Worker completed - marker_idx: {marker_idx}, successful: {successful_trees}/{len(trees_in_marker)}, "
              f"cpu_core: {cpu_core}, node: {node_name}, "
              f"start_time: {start_time}, end_time: {end_time}, "
              f"duration: {marker_duration} minutes", flush=True)'''
        
        return results



        
    def reconstruct_single_gene_tree_for_subtree_for_astral(self, tree_id, data, alignment_idx, alignment):
        """
        Reconstruct a single ML tree for one gene for one subtree.
        Returns the tree_id, alignment_idx, and the tree string.
        """
        start = time.time()
        cpu_core_for_each_subtree = self.get_cpu_core()
        node_name = platform.node()
        start_time = datetime.datetime.now()
        
        #print(f'cpu_core: {cpu_core}, node: {node_name}', flush=True)
        
        # Create node and extract alignment
        node = MyNode(
            data=data["bitmap"],
            data_with_outgroup=data["bitmap_with_outgroup"],
            outgroup_name=data["outgroup_name"],
            tree_id=tree_id
        )
        
        # Select rows for this subtree
        selected_leaves_names = [leave_name for leave_name, flag in zip(self.leaves_names, node.bitmap_with_outgroup) if flag]
        selected_taxa_metadata = spectraltree.TaxaMetadata(dendropy.TaxonNamespace(selected_leaves_names), selected_leaves_names)
        alignment_selected_rows = alignment
        
        # Reconstruct ML tree
        raxml_instant = spectraltree.RAxML()
        sub_ml_tree = raxml_instant(alignment_selected_rows, selected_taxa_metadata)
        tree_string = sub_ml_tree.as_string('newick')
        
        # Timing info
        end_time = time.time()
        duration = round((end_time - start) / 60, 2)
        print(f'- tree_id: {tree_id}, alignment_idx: {alignment_idx}, cpu_core: {cpu_core_for_each_subtree}, pid: {os.getpid()}, node: {node_name}, start_time: {start_time}, end_time: {datetime.datetime.now()}, completed in {duration} minutes', flush=True)
        print()
            
        return tree_id, alignment_idx, tree_string

    def reconstruct_small_subtrees_parallel(self):

        folder_export_subtrees = self.create_folder_to_subtrees()

        """
        Reconstruct trees in parallel for nodes with small enough taxon count.
        Stores results back in self.nodes_data. - can be done only on the cluster.
        """
        self.parallel_tree_infos = []


        #estimated_peak_memory_for_trees = self.peak_memory_similarity * 1.52
        estimated_peak_memory_for_trees_in_gb = self.number_of_genes / 7
        #print('estimated_peak_memory_for_trees_in_gb acording number of genes')
        #print(estimated_peak_memory_for_trees_in_gb)
        estimated_cpu_number_for_trees = np.floor(1000 / estimated_peak_memory_for_trees_in_gb).astype(int)
        print('estimated_cpu_number_for_trees')
        print(estimated_cpu_number_for_trees)
        cpus_on_node = int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 1))
        #self.max_workers = min(estimated_cpu_number_for_trees, cpus_on_node) 

        self.max_workers = int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 1))
        print(f"Using {self.max_workers} CPU workers for tree reconstruction")
    
        # Filter only the nodes that should be reconstructed
        self.nodes_to_reconstruct = {
            tree_id: data
            for tree_id, data in self.nodes_data.items()
            if np.sum(data["bitmap"]) <= self.threshold
        }
        #print('self.nodes_to_reconstruct befor:')
        #print(self.nodes_to_reconstruct)

        # Sort items by descending bitmap sum
        self.nodes_to_reconstruct = dict(
            sorted(
                self.nodes_to_reconstruct.items(),
                key=lambda item: np.sum(item[1]["bitmap"]),
                reverse=True
            )
        )
        
        print(f'number of trees need recunstruction: {len(self.nodes_to_reconstruct)}')
        
        #print(f'******peak_memory in kb befor all recunstructions {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')
        
        #peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        #peak_gb = peak_kb / 1024 / 1024
        #print(f'***peak_memory GB befor reconstruction :{peak_gb}***')
        
        if self.inner_method == 'astral':

            folder_export_subtrees = self.create_folder_to_subtrees()
            
            # Handle different number of genes for reconstruction
            if different_num_of_genes_to_recunstruct:
                self.alignments = self.alignments[:number_of_genes_for_reconstruction]
                print(f'Partition run on {self.number_of_genes} and reconstruction run on {len(self.alignments)} genes')
        
            # Phase 1: Reconstruct all gene trees in parallel (in chunks)
            print("=== Phase 1: Reconstructing all gene trees in parallel ===", flush=True)

            # Dictionary to collect gene trees: {tree_id: [tree_strings]}
            gene_trees_by_subtree = {tree_id: [] for tree_id in self.nodes_to_reconstruct.keys()}

            # Create list of all tasks
            all_tasks = []
            for alignment_idx, alignment in enumerate(self.alignments):
                alignment_list = []
                for tree_id, data in self.nodes_to_reconstruct.items():
                    # Extract the alignment rows for this subtree
                    bitmap_with_outgroup = data["bitmap_with_outgroup"]
                    alignment_selected_rows = alignment[bitmap_with_outgroup, :]
                    
                    # Append as tuple: (tree_id, data, alignment_idx, alignment_subset)
                    alignment_list.append((tree_id, data, alignment_idx, alignment_selected_rows))
                all_tasks.append(alignment_list)
            
            total_tasks = len(all_tasks)
            print(f"Total tasks to process: {total_tasks}", flush=True)
            
            # Calculate chunk size: total tasks / max workers
            chunk_size = len(self.nodes_to_reconstruct)  # Number of subtrees per chunk
            num_of_marker = len(self.alignments)            # Number of genes/chunks
            
            print(f"Chunk size: {chunk_size}, Number of chunks: {num_of_marker}", flush=True)
            
            # Process chunks in parallel (each chunk = one gene with all its subtrees)
            gene_trees_per_subtree = {tree_id: [] for tree_id in self.nodes_to_reconstruct.keys()}
            
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                # Submit each chunk to a separate CPU
                # Each chunk contains all subtrees for one gene
                chunk_futures = []
                
                for marker_idx, trees_in_marker in enumerate(all_tasks):
                    print(f"Submitting chunk {marker_idx + 1}/{num_of_marker} ({len(trees_in_marker)} subtrees)", flush=True)
                    
                    # Submit the entire chunk to one worker
                    future = executor.submit(
                        self.reconstruct_all_trees_for_single_marker_for_astral,
                        marker_idx,
                        trees_in_marker
                    )
                    chunk_futures.append((marker_idx, future))
                
                # Collect results as chunks complete
                for marker_idx, future in chunk_futures:
                    chunk_results = future.result()
                    
                    # Merge results: chunk_results = {tree_id: tree_string, ...}
                    for tree_id, tree_string in chunk_results.items():
                        if tree_string is not None:
                            gene_trees_per_subtree[tree_id].append(tree_string)
                    
                    print(f"Chunk {marker_idx + 1}/{num_of_marker} completed", flush=True)
                        
            
            print(f"All chunks processed!", flush=True)
    

            # write to the disk for astral
            for tree_id, gene_trees_list in gene_trees_per_subtree.items():
                # Write all gene trees for this subtree to a file
                ml_trees_file = os.path.join(folder_export_subtrees, f"{tree_id}_ml_trees.newick")
                with open(ml_trees_file, 'w') as f:
                    for tree_string in gene_trees_list:
                        f.write(tree_string)
            

            # Phase 2: Run ASTRAL in parallel for each subtree
            print("=== Phase 2: Running ASTRAL in parallel for each subtree ===", flush=True)
            
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                futures = []
                
                for tree_id, data in self.nodes_to_reconstruct.items():
                    future = executor.submit(
                        self.run_astral_for_subtree,
                        tree_id, data, folder_export_subtrees
                    )
                    futures.append(future)
                
                # Collect results and update nodes_data
                for future in as_completed(futures):
                    tree_id, sub_tree = future.result()
                    self.nodes_data[tree_id]["tree"] = sub_tree
                    #print(f'Tree saved to nodes_data for {tree_id}', flush=True)
        

        else:
        
        
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                for tree_id, data in self.nodes_to_reconstruct.items():
                    executor.submit(self._reconstruct_single_tree, tree_id, data)
            
            # retrive back all the trees:
            for tree_id, data in self.nodes_data.items():
                if np.sum(data["bitmap"]) <= self.threshold:
                    subtree_file_for_parrelal = os.path.join(folder_export_subtrees, f"{tree_id}.nwk")
                    while not os.path.exists(subtree_file_for_parrelal):
                        time.sleep(0.5)  # Sleep half a second before checking again
                    # File now exists; load it
                    with open(subtree_file_for_parrelal, 'r') as f:
                        loaded_subtree = Tree.get(file=f, schema="newick")
                    self.nodes_data[tree_id]["tree"] = loaded_subtree
                
        #print('final nodes that recunstructed data after:')
        #print(self.nodes_data)



        # Delete the directory and its contents
        shutil.rmtree(folder_export_subtrees)

        
        #print(f'******peak_memory in kb after all recunstructions {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}')



        
    

    

    # call all functions
    def deep_spectral_tree_reconstruction(self, verbose = False):
        total_time_pre_processing = 0
        total_time_create_similarity = 0
        total_time_merging = 0
        total_time_partitioning = 0
        total_time_recunstructing = 0


        start_time_pre_processing = time.time()
        
        self.pre_process()

        end_time_pre_processing = time.time()
        time_pre_processing = (end_time_pre_processing - start_time_pre_processing) / 60
        total_time_pre_processing += time_pre_processing
        

        start_time_create_similarity = time.time()

        self.create_similarity_matrices()
        
        end_time_create_similarity = time.time()
        time_create_similarity = (end_time_create_similarity - start_time_create_similarity) / 60
        total_time_create_similarity += time_create_similarity
        

        
        start_time_partitioning = time.time()

        partitioning_tree = MyTree([True]*len(self.taxa_metadata), [True]*len(self.taxa_metadata))
        cur_node = partitioning_tree.root
        
        make_partition = 'yes'
        self.make_partition = make_partition

        self.nodes_data = {}  # key = tree_id, value = dict with bitmap info and outgroup
        
        while True:

            if (cur_node.right != None) and (cur_node.right.tree !=None) and (cur_node.left.tree != None):
                
                cur_node.tree = 'yes'

                if cur_node.parent == None:
                    break
                if cur_node.parent.right == cur_node:
                    cur_node = cur_node.parent.left
                else:
                    cur_node = cur_node.parent

                
            elif np.sum(cur_node.bitmap) > self.threshold:
                
                bitmap1, bitmap2, bitmap_with_outgroup_1, bitmap_with_outgroup_2, outgroup_name_1, outgroup_name_2 = self.splitTaxa(cur_node)#,num_gaps,min_split)
                if verbose: print("partition")
                if verbose: print("L1 size: ", np.sum(bitmap1))
                if verbose: print("L2 size: ", np.sum(bitmap2))
                cur_node.setLeft(MyNode(bitmap1, bitmap_with_outgroup_1, outgroup_name_1))
                cur_node.setRight(MyNode(bitmap2, bitmap_with_outgroup_2, outgroup_name_2))

                # Save info for each
                self.nodes_data[cur_node.left.tree_id] = {
                    "bitmap": cur_node.left.bitmap,
                    "bitmap_with_outgroup": cur_node.left.bitmap_with_outgroup,
                    "outgroup_name": cur_node.left.outgroup_name
                }
                self.nodes_data[cur_node.right.tree_id] = {
                    "bitmap": cur_node.right.bitmap,
                    "bitmap_with_outgroup": cur_node.right.bitmap_with_outgroup,
                    "outgroup_name": cur_node.right.outgroup_name
                }

                cur_node = cur_node.right

            else:
                #print(cur_node.tree_id)
                
                cur_node.tree = 'yes'

                if cur_node.parent == None:
                    break
                if cur_node.parent.right == cur_node:
                    cur_node = cur_node.parent.left
                else:
                    cur_node = cur_node.parent

        
        # add data for the root
        self.nodes_data[''] = {
            "bitmap": [a or b for a, b in zip(self.nodes_data['l']['bitmap'], self.nodes_data['r']['bitmap'])],
            "bitmap_with_outgroup": [a or b for a, b in zip(self.nodes_data['l']['bitmap_with_outgroup'], self.nodes_data['r']['bitmap_with_outgroup'])],
            "outgroup_name": None}

        end_time_partitioning = time.time()
        time_partitioning = (end_time_partitioning - start_time_partitioning) / 60
        total_time_partitioning += time_partitioning

        
        del self.similarities_array
        del self.partial_fiedler
        del self.partial_partition
        

        if self.type_of_run == 'cluster':
            # find the memory usage i use
            process = psutil.Process(os.getpid())
            current_memory_gb = process.memory_info().rss / 1024 / 1024 / 1024
            print(f'***Current memory GB after similarity and befor reconstruction: {current_memory_gb:.2f}***')

        start_time_recunstructing = time.time()
        
        # Reconstruct trees for all nodes in nodes_data with few enough taxa
        if self.type_of_run == 'pc':
            self.reconstruct_small_subtrees()
        if self.type_of_run == 'cluster':            
            self.reconstruct_small_subtrees_parallel()


        end_time_recunstructing = time.time()
        time_recunstructing = (end_time_recunstructing - start_time_recunstructing) / 60
        total_time_recunstructing += time_recunstructing

        
        start_time_merging = time.time()

        estimated_tree = self.merge_trees_for_parallel(self.nodes_data)

        end_time_merging = time.time()
        time_merging = (end_time_merging - start_time_merging) / 60
        total_time_merging += time_merging


        print(estimated_tree.as_ascii_plot())

        self.export_tree(estimated_tree)

        return estimated_tree, total_time_pre_processing, total_time_create_similarity, total_time_partitioning, total_time_merging, total_time_recunstructing


    # plain reconstruction (without partitions)
    def tree_reconstruction(self):
        self.pre_process()

        partitioning_tree = MyTree([True]*len(self.taxa_metadata), [True]*len(self.taxa_metadata))
        cur_node = partitioning_tree.root

        # reminder to the code after not to do re root with outgroup where there is no 2 subtrees
        make_partition = 'no'
        self.make_partition = make_partition
        
        cur_node.tree = self.reconstruct_alg_wrapper(cur_node)#, **kargs)

        estimated_tree = partitioning_tree.root.tree
        self.export_tree(estimated_tree)
        return estimated_tree

    def create_true_partitions(self, true_subtree):
    
        # Export the pruned tree to a Newick string
        #true_tree_newick_str = true_tree.as_string(schema="newick")
        # Read the Newick string into a Biopython Phylo tree
        #true_tree = Phylo.read(StringIO(true_tree_newick_str), "newick")
    
        root = true_subtree.root
    
        
        #find the name of all the leaves
        leaves = true_subtree.get_terminals()
    
        leaves_names =[]
        # Print the names of the leaf nodes
        for leaf in leaves:
            name = leaf.name
            leaves_names.append(name)
        leaves_names.sort()
        #print('leaves_names.sort()')
        #print(leaves_names)
              
    
        #divied to 2 groups
        branch_true_1, branch_true_2 = root.clades
        sub_tree_1 = Phylo.BaseTree.Tree(branch_true_1)
        sub_tree_2 = Phylo.BaseTree.Tree(branch_true_2)
        
        #find the name of sub group 1
        leaves_sub_tree_1 = sub_tree_1.get_terminals()
    
        leaves_names_subtree_1 =[]
        # Print the names of the leaf nodes
        for leaf in leaves_sub_tree_1:
            name = leaf.name
            leaves_names_subtree_1.append(name)
        leaves_names_subtree_1.sort()
    
        #find the name of sub group 2
        leaves_sub_tree_2 = sub_tree_2.get_terminals()
    
        leaves_names_subtree_2 =[]
        # Print the names of the leaf nodes
        for leaf in leaves_sub_tree_2:
            name = leaf.name
            leaves_names_subtree_2.append(name)
        leaves_names_subtree_2.sort()
        
        #create true filder according true tree and leaves list
    
        true_fiedler = []
        for name in leaves_names:
            if name in leaves_names_subtree_1:
                true_fiedler.append(1)
            else:
                true_fiedler.append(0)
            
            
        # Find all clades in the tree
        clades = list(root.find_clades())
    
        # Get the terminal nodes of each clade
        clades_groups = []
        for clade in clades:
            terminals = clade.get_terminals()
            clade_names = []
            for terminal in terminals:
                name = terminal.name
                clade_names.append(name)
            clades_groups.append(clade_names)
    
        #ommit all the groups with 1 terminal node
        clades_groups_without_one_terminal_nodes = []
        for i in list(range(0,len(clades_groups))):
            group = clades_groups[i]
            if ((len(group)>1) and (len(group)<len(leaves_names))):
                clades_groups_without_one_terminal_nodes.append(group)
    
        # form true fidler vector to each clade
        all_true_fiedlers = []
        for i in list(range(0,len(clades_groups_without_one_terminal_nodes))):
            each_true_fiedler = []
            for name in leaves_names:
                if name in clades_groups_without_one_terminal_nodes[i]:
                    each_true_fiedler.append(1)
                else:
                    each_true_fiedler.append(0)
            all_true_fiedlers.append(each_true_fiedler)
    

        # Create a new list with the flipped sublists (0 -> 1 and 1 -> 0)
        flipped_fiedlers = [[1 - x for x in sublist] for sublist in all_true_fiedlers]
        
        # Add the flipped sublists to the original list
        all_true_fiedlers.extend(flipped_fiedlers)
        
        #self.all_true_fiedlers = all_true_fiedlers
        #print(self.all_true_fiedlers)
        return all_true_fiedlers


    def find_score(self, all_true_fiedlers, estimated_fiedler):
        #print('all_true_fiedlers')
        #print(all_true_fiedlers)
        #print('estimated_fiedler')
        #print(estimated_fiedler)
        #find all the scores
        all_score = []
        for i in list(range(0,len(all_true_fiedlers))):
            '''
            if (np.count_nonzero(estimated_fiedler == 1)<=1) or (np.count_nonzero(estimated_fiedler == 0)<=1) or (np.count_nonzero(estimated_fiedler == 1)>=len(estimated_fiedler)-1) or (np.count_nonzero(estimated_fiedler == 0)>=len(estimated_fiedler)-1):
                score = 0
            else:'''
            score = rand_score(all_true_fiedlers[i], estimated_fiedler)
            #print('score')
            #print(score)
            all_score.append(score)
            #print('all_score')
            #print(all_score)
        final_score = max(all_score)
        
        return final_score

    def save_only_higher_value_from_each_pair(self, data):
        # save only the higher value from each pair
        result = {}
        checked = set()
        
        for key in data:
            # Only consider keys that are not the root (empty string)
            if key == '':
                continue
            parent = key[:-1]
            if parent + 'l' in data and parent + 'r' in data:
                left_key = parent + 'l'
                right_key = parent + 'r'
                if (left_key, right_key) in checked or (right_key, left_key) in checked:
                    continue  # Already processed this pair
                # Keep the key with the higher value
                if data[left_key] >= data[right_key]:
                    result[left_key] = data[left_key]
                else:
                    result[right_key] = data[right_key]
                checked.add((left_key, right_key))
        
        result.pop('', None)
        
        return result
        

    def add_parent_true_leaves(self):
        for node_key, node_data in self.nodes_data.items():
            # Determine parent key
            if node_key == '':
                parent_key = None
            else:
                parent_key = node_key[:-1]
            
            # Get parent's bitmap if exists
            if parent_key in self.nodes_data:
                parent_bitmap = self.nodes_data[parent_key]['bitmap']
                # Get leaves corresponding to True in parent's bitmap
                parent_true_leaves = [
                    leaf for leaf, bit in zip(self.leaves_names, parent_bitmap) if bit
                ]
            else:
                parent_true_leaves = []
                
            if node_key == '':
                parent_true_leaves = []

            # Add to node data
            node_data['parent_true_leaves'] = sorted(parent_true_leaves)

    def add_leaves(self):
        for node_key, node_data in self.nodes_data.items():
            
        # Get parent's bitmap if exists
            node_bitmap = self.nodes_data[node_key]['bitmap']
            # Get leaves corresponding to True in parent's bitmap
            node_leaves = [
                leaf for leaf, bit in zip(self.leaves_names, node_bitmap) if bit
            ]

            # Add to node data
            node_data['node_leaves'] = node_leaves
            

    
    def cut_subtrees_by_parent_leaves(self):
        for node_key, node_data in self.nodes_data.items():
            parent_leaves = node_data.get('parent_true_leaves', [])
            if not parent_leaves:
                continue  # Nothing to cut, skip
            
            # Deep copy the tree to avoid modifying the original
            sub_tree = copy.deepcopy(true_tree)
            
            # Prune all leaves NOT in parent_leaves
            to_prune = [term for term in sub_tree.get_terminals() if term.name not in parent_leaves]
            for term in to_prune:
                sub_tree.prune(term)

            # Store the subtree (could also export as Newick, etc.)
            node_data['parent_subtree'] = sub_tree
            #print(node_data['parent_subtree'])

            
    def find_true_partition(self, verbose = False):
        #import true tree and remove outgroup
        #true_tree = dendropy.Tree.get(file=open(true_tree_path, 'r'), schema="newick")
        #true_tree.prune_taxa_with_labels([outgroup_name])
        
        #find all the true partitions in the true tree
        #self.create_true_partitions()
        #print(self.all_true_fiedlers)


        total_time_pre_processing = 0
        total_time_create_similarity = 0
        total_time_merging = 0
        total_time_partitioning = 0
        total_time_recunstructing = 0

        start_time_pre_processing = time.time()
        
        self.pre_process()

        end_time_pre_processing = time.time()
        time_pre_processing = (end_time_pre_processing - start_time_pre_processing) / 60
        total_time_pre_processing += time_pre_processing
        
        start_time_create_similarity = time.time()

        self.create_similarity_matrices()
        
        end_time_create_similarity = time.time()
        time_create_similarity = (end_time_create_similarity - start_time_create_similarity) / 60
        total_time_create_similarity += time_create_similarity

        
        start_time_partitioning = time.time()

        partitioning_tree = MyTreeForPartitionCheck([True]*len(self.taxa_metadata))
        cur_node = partitioning_tree.root
        
        make_partition = 'yes'
        self.make_partition = make_partition

        self.nodes_data = {}  # key = tree_id, value = dict with bitmap info and outgroup
        
        while True:

            if (cur_node.right != None) and (cur_node.right.tree !=None) and (cur_node.left.tree != None):
                
                cur_node.tree = 'yes'

                if cur_node.parent == None:
                    break
                if cur_node.parent.right == cur_node:
                    cur_node = cur_node.parent.left
                else:
                    cur_node = cur_node.parent

                
            elif np.sum(cur_node.bitmap) > self.threshold:
                
                bitmap1, bitmap2 = self.splitTaxaForPartitionCheck(cur_node)#,num_gaps,min_split)
                if verbose: print("partition")
                if verbose: print("L1 size: ", np.sum(bitmap1))
                if verbose: print("L2 size: ", np.sum(bitmap2))
                cur_node.setLeft(MyNodeForPartitionCheck(bitmap1))
                cur_node.setRight(MyNodeForPartitionCheck(bitmap2))

                # Save info for each
                self.nodes_data[cur_node.left.tree_id] = {
                    "bitmap": cur_node.left.bitmap,
                }
                self.nodes_data[cur_node.right.tree_id] = {
                    "bitmap": cur_node.right.bitmap,
                }

                cur_node = cur_node.right
                #print(self.leaves_names)

            else:
                #print(cur_node.tree_id)
                
                cur_node.tree = 'yes'

                if cur_node.parent == None:
                    break
                if cur_node.parent.right == cur_node:
                    cur_node = cur_node.parent.left
                else:
                    cur_node = cur_node.parent

        
        # add data for the root
        self.nodes_data[''] = {
            "bitmap": [a or b for a, b in zip(self.nodes_data['l']['bitmap'], self.nodes_data['r']['bitmap'])]}

        
        self.add_parent_true_leaves()
        self.add_leaves()
        #print('**self.nodes_data')
        #print(self.nodes_data)
        #print('***')
        
        #print('self.nodes_data:')
        #print(self.nodes_data)
        '''
        # 1. Filter only the nodes that should be reconstructed
        self.nodes_data = {
            tree_id: data
            for tree_id, data in self.nodes_data.items()
            if np.sum(data["bitmap"]) <= self.threshold
        }'''
        
        # 2. For each item, cut a sub tree according to parent leaves
        for node_key, node_data in self.nodes_data.items():
            parent_leaves = node_data['parent_true_leaves']

            
            # Deep copy the tree to avoid modifying the original
            sub_tree = copy.deepcopy(true_tree)
        
            # Prune all leaves NOT in parent_leaves
            #if parent_leaves is not None:
            if node_key == '':
                sub_tree = ''
            else:
                to_prune = [term for term in sub_tree.get_terminals() if term.name not in parent_leaves]
                for term in to_prune:
                    sub_tree.prune(term)
        
            # Store the subtree (could also export as Newick, etc.)
            node_data['parent_subtree'] = sub_tree
        
        # 3. Now, create new dict for filtered nodes, copying all fields
        node_data = {}
        
        for tree_id, content in self.nodes_data.items():
            if tree_id == '':
                continue
            bitmap_10 = [1 if x else 0 for x in content['bitmap']]
            node_data[tree_id] = {
                'bitmap_10': bitmap_10,
                'ri': None,
                'bitmap': content['bitmap'],
                'node_leaves': content['node_leaves'],
                'parent_true_leaves': content['parent_true_leaves'],
                'parent_subtree': content.get('parent_subtree', None),   # <-- Now this will be the subtree!
            }
        
        self.nodes_data = node_data
        

        for tree_id, content in self.nodes_data.items():
            subtree = content['parent_subtree']
            # If you are certain all valid subtrees are not None, you don't need the .get()
            all_true_fiedlers = self.create_true_partitions(subtree)
            content['all_true_fiedlers'] = all_true_fiedlers


        for tree_id, content in self.nodes_data.items():
            node_leaves = content['node_leaves']
            parent_true_leaves = content['parent_true_leaves']
            # Build bitmap: 1 if node_leaves contains the parent leaf, else 0
            node_bitmap_in_parent = [1 if leaf in node_leaves else 0 for leaf in parent_true_leaves]
            content['estimated_fiedler'] = node_bitmap_in_parent


        for tree_id, content in self.nodes_data.items():
            all_true_fiedlers = content.get('all_true_fiedlers')
            estimated_fiedler = content.get('estimated_fiedler')

            # Compute the score
            ri_score = self.find_score(all_true_fiedlers, estimated_fiedler)
            content['ri_score'] = ri_score
        '''
        print('node_data_after_building:')
        for tree_id, content in self.nodes_data.items():
                print(f"tree_id: {tree_id}")
                print(f"node_leaves: {content['node_leaves']}")
                print(f"node_leaves_len: {len(content['node_leaves'])}")
                print(f"parent_true_leaves: {content['parent_true_leaves']}")
                print(f"parent_true_leaves_len: {len(content['parent_true_leaves'])}")

                print("------")'''
        #print(self.nodes_data)

        summary = []
        for tree_id, node in self.nodes_data.items():
            summary.append({
                'tree_id': tree_id,
                'bitmap_sum': sum(node['bitmap']),
                'ri_score': node.get('ri_score', None)  # use .get in case not present
            })
        
        #print(summary)
        #print('****')
        return summary
        



        #self.ri_scores = sum(self.ri_scores.values()) / len(self.ri_scores)
        
        end_time_partitioning = time.time()
        time_partitioning = (end_time_partitioning - start_time_partitioning) / 60
        total_time_partitioning += time_partitioning        

        #return self.ri_scores



class MyNode(object):
    def __init__(self, data, data_with_outgroup, outgroup_name, tree_id=''):
        self.bitmap = data
        self.bitmap_with_outgroup = data_with_outgroup
        self.outgroup_name = outgroup_name
        self.tree = None
        self.left = None
        self.right = None
        self.parent = None
        self.tree_id = tree_id

    def setLeft(self,node):
        self.left = node
        node.parent = self
        node.tree_id = self.tree_id + 'l'  # Add 'l' to path
        
    def setRight(self,node):
        self.right = node
        node.parent = self
        node.tree_id = self.tree_id + 'r'  # Add 'r' to path
        
class MyTree(object):
    def __init__(self, data = None, data_with_outgroup = None, outgroup_name = None):
        self.root = MyNode(data, data_with_outgroup, outgroup_name, tree_id='')  # Root has empty path
    def setLeft(self,node):
        self.left = node
        node.parent = self
    def setRight(self,node):
        self.right = node
        node.parent = self


class MyNodeForPartitionCheck(object):
    def __init__(self, data, tree_id=''):
        self.bitmap = data
        self.tree = None
        self.left = None
        self.right = None
        self.parent = None
        self.tree_id = tree_id

    def setLeft(self,node):
        self.left = node
        node.parent = self
        node.tree_id = self.tree_id + 'l'  # Add 'l' to path
        
    def setRight(self,node):
        self.right = node
        node.parent = self
        node.tree_id = self.tree_id + 'r'  # Add 'r' to path
        
class MyTreeForPartitionCheck(object):
    def __init__(self, data = None):
        self.root = MyNodeForPartitionCheck(data, tree_id='')  # Root has empty path
    def setLeft(self,node):
        self.left = node
        node.parent = self
    def setRight(self,node):
        self.right = node
        node.parent = self

