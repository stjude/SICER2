# Author: 2022 STJUDE Modupeore Adetunji

# Separate bed file to individual chromosomes

import os
import re
import sys
import subprocess
from functools import partial
import numpy as np

def separate_bedpe_chroms(file, chrom):
    file_name = os.path.basename(file)
    file_name = file_name.replace('.bed', '')

    print_return = ""

    match = chrom + "[[:space:]]"
    individual_bed = subprocess.Popen(['grep', match, file], stdout=subprocess.PIPE)
    bed_reads = str(individual_bed.communicate()[0], 'utf-8').splitlines()
    read_dtype = np.dtype([('chrom', 'U6'), ('start', np.int32), ('end', np.int32), ('name', 'U20'), ('score', np.int32), ('strand', 'U1')])
    processed_reads = np.empty(len(bed_reads), dtype=read_dtype)

    reads_count = 0
    for i, reads in enumerate(bed_reads):
        reads = re.split('\t', reads)
        if len(reads) < 6:
            sys.stderr.write(
                "Error: Input BED files must have the first six fields. Check " + file_name + " to see if it has the following fields: chrom, chromStart, chromEnd, name, score, and strand\n")
            sys.exit(1)
        reads[1] = int(reads[1])
        reads[2] = int(reads[2])
        processed_reads[i] = tuple(reads)
        reads_count += 1

    print_return += ('{:<5s}{:^25d}'.format(chrom, reads_count))
    name_for_save = file_name + "_" + chrom + ".npy"
    np.save(name_for_save, processed_reads)

    return (print_return, reads_count)

def main(args, file, pool):
    chroms = args.species_chroms

    separate_chroms_partial = partial(separate_bedpe_chroms, file)
    results_count = pool.map(separate_chroms_partial, chroms)

    total_read_count = 0
    print(('-' *30))
    print(('{:<5s}{:^25s}'.format("chrom", "Total reads")))
    print(('-' *30))
    for result in results_count:
        print(result[0])
        total_read_count += result[1]

    return total_read_count
