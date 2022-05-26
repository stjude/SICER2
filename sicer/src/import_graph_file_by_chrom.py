# Written: 2022 STJUDE Modupeore Adetunji

from functools import partial
import os
import re
import subprocess
import sys
import numpy as np

def match_chrom(file, chrom):
    match = chrom + "[[:space:]]"
    matched_reads = subprocess.Popen(['grep', match, file], stdout=subprocess.PIPE) #Use Popen so that if no matches are found, it doesn't throw an exception
    chrom_reads = str(matched_reads.communicate()[0], 'utf-8').splitlines()  # generates a list of each reads, which are represented by a string value
    file_name = os.path.basename(file)
    chrom_graph = []
    chrom_data = []
    tag_count = 0
    for i, reads in enumerate(chrom_reads):
        reads = re.split('\t', reads)
        if len(reads) < 4:
            sys.stderr.write(
                "Error: Input files must have the first six fields. Check " + file_name + " to see if it has the following fields: chrom, chromStart, chromEnd, and score\n")
            sys.exit(1)
        reads[1] = int(reads[1])
        reads[2] = int(reads[2])
        reads[3] = int(reads[3])
        if reads[3] > 0:
            output1 = reads
            output2 = (chrom, reads[1], reads[2], reads[3])
            chrom_data.append(output1)
            chrom_graph.append(output2)
            tag_count += reads[3]

    return (chrom_data, chrom_graph, tag_count)


def separate_reads(path_to_file, chrom):
    file_name = os.path.basename(path_to_file)
    file_name = file_name.replace('.bed', '')

    (chrom_reads, chrom_graph, tag_count) = match_chrom(path_to_file, chrom)  # Separates all reads by chromosome

    file_save_name = file_name + '_' + chrom + '.npy'
    graph_save_name = file_name + '_' + chrom + '_graph.npy'

    np_chrom_graph = np.array(chrom_graph, dtype=object)
    np.save(graph_save_name, np_chrom_graph)

    np_chrom_reads = np.array(chrom_reads, dtype=object)
    np.save(file_save_name, np_chrom_reads)

    return tag_count


def main(args, path_to_file, pool):
    chroms = args.species_chroms

    separate_reads_partial = partial(separate_reads, path_to_file)
    tag_counts = pool.map(separate_reads_partial, chroms)

    total_tag_count = 0
    for result in tag_counts:
        total_tag_count += result

    return total_tag_count
