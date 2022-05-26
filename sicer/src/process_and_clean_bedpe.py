# Written: 2022 STJUDE Modupeore Adetunji

from functools import partial
import subprocess
import re
import numpy as np

def graph_bins_chrom(file, chrom):
    file_name = file.replace('.bed', '')
    bed_file_name = file_name + '_' + chrom + '.npy'
    new_file_name = file_name + '_' + chrom + '-results.bed'

    print_return = ""
    chrom_reads = np.load(bed_file_name, allow_pickle=True)
    with open(new_file_name, 'w') as outfile:
        for read in chrom_reads:
            chrom = read[0]
            start = read[1]
            end = read[2]
            name = read[3]
            score = read[4]
            strand = read[5]
            line = (chrom + '\t' + str(start) + '\t' + str(end) + '\t' + name + '\t' + str(score) + '\t' + strand + '\n')
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
            output1 = (chrom, reads[1], reads[2], 'placeholder', '255', '+')
            output2 = (chrom, reads[1], reads[2], reads[3])
            chrom_data.append(output1)
            chrom_graph.append(output2)
            tag_count += reads[3]

    file_save_name = file_name + '_' + chrom + '.npy'
    graph_save_name = file_name + '_' + chrom + '_graph.npy'

    np_chrom_graph = np.array(chrom_graph, dtype=object)
    np.save(graph_save_name, np_chrom_graph)

    np_chrom_reads = np.array(chrom_data, dtype=object)
    np.save(file_save_name, np_chrom_reads)

    print_return += ('Total count of ' + chrom + ' tags: ' + str(tag_count) + '\n')
    return (print_return, tag_count)

def main(args, file, pool):
    chroms = args.species_chroms

    graph_bins_chrom_partial = partial(graph_bins_chrom, file)
    tag_counts = pool.map(graph_bins_chrom_partial, chroms)
    total_tag_count = 0
    print_return = ""
    for result in tag_counts:
        print_return += result[0]
        #print(result[0])
        total_tag_count += result[1]

    return (print_return, total_tag_count)
