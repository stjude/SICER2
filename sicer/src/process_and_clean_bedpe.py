from functools import partial
import subprocess
import os
import sys
import re
import numpy as np

def graph_bins_chrom(file, chrom):
    file_name = os.path.basename(file)
    file_name = file_name.replace('.bed', '')
    bed_file_name = file_name + '_' + chrom + '.npy'
    new_file_name = file_name + '_' + chrom + '-results.bed'

    chrom_reads = np.load(bed_file_name, allow_pickle=True)
    with open(new_file_name, 'w') as outfile:
        for read in chrom_reads:
            chrom = read[0]
            start = read[1]
            end = read[2]
            name = read[3]
            score = read[4]
            strand = read[5]
            line=(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + name + '\t' + str(score) + '\t' + strand + '\n')
            outfile.write(line)

    windows = chrom + '.windows'
    graph_reads = subprocess.Popen(['bedtools', 'intersect', '-c', '-a', windows, '-b', new_file_name], stdout=subprocess.PIPE)
    chrom_reads = str(graph_reads.communicate()[0], 'utf-8').splitlines()

    chrom_graph = []
    chrom_data = []
    tag_count = 0
    for i, reads in enumerate(chrom_reads):
        reads = re.split('\t', reads)
        reads[1] = int(reads[1])
        reads[2] = int(reads[2])
        reads[3] = int(reads[3])
        if reads[3] > 0:
            output1 = reads
            output2 = (chrom, reads[1], reads[2], reads[3])
            chrom_data.append(output1)
            chrom_graph.append(output2)
            tag_count += reads[3]

    file_save_name = file_name + '_' + chrom + '.npy'
    graph_save_name = file_name + '_' + chrom + '_graph.npy'

    np_chrom_graph = np.array(chrom_graph, dtype=object)
    np.save(graph_save_name, np_chrom_graph)

    np_chrom_reads = np.array(chrom_reads, dtype=object)
    np.save(file_save_name, np_chrom_reads)

    print('Total count of ' + chrom + ' tags: ' + str(tag_count))
    return tag_count

def main(args, file, pool):
    chroms = args.species_chroms

    graph_bins_chrom_partial = partial(graph_bins_chrom,file)
    tag_counts = pool.map(graph_bins_chrom_partial, chroms)
    total_tag_count = 0
    for result in tag_counts:
        total_tag_count += result

    return total_tag_count
