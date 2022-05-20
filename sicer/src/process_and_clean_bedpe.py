from functools import partial
import subprocess
import os
import sys
import re
import numpy as np

def separate_bedpe_chroms(file, chrom):
    file_name = os.path.basename(file)
    file_name = file_name.replace('.bedpe', '')

    match = chrom + "[[:space:]]"
    syntax = 'grep %s %s > %s_%s-results.bed' % (match, file, file_name, chrom)
    subprocess.Popen(syntax, stdin=subprocess.PIPE, shell=True,)

    print("Cleaning BEDPE file for chrom %s..." % (chrom))
    awkcommand = 'awk -F\\\\t \'{if($1==$4 && $1~"chr" && $6-$2<1000) { print $1 "\\t" $2 "\\t" $6 "\\t" $7 } }\''

    new_file_name = file_name + '_' + chrom
    process = subprocess.Popen("%s %s-results.bed > %s-cleaned_bedpe.bed" % (awkcommand, new_file_name, new_file_name), stdin=subprocess.PIPE, shell=True,)
    process.communicate()
    if process.returncode != 0:
        sys.stderr.write("Error: Cannot clean bedpe file.\nCheck if bedtools2 (https://github.com/arq5x/bedtools2) has been installed correctly.\n")
        sys.exit(1)

    windows = chrom + '.windows'
    cleaned = new_file_name + '-cleaned_bedpe.bed'
    graph_reads = subprocess.Popen(['intersectBed', '-c', '-a', windows, '-b', cleaned], stdout=subprocess.PIPE)
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

    return tag_count


def main(args, file, pool):
    chroms = args.species_chroms
    separate_chroms_partial = partial(separate_bedpe_chroms, file)

    tag_counts = pool.map(separate_chroms_partial, chroms)
    for result in tag_counts:
        total_tag_count += result

    return total_tag_count
