# Written: 2022 STJUDE Modupeore Adetunji

from functools import partial
import subprocess
import sys

def create_bed_windows(chrom, chrom_length, bin_size):
    file_save_name = chrom + '.windows'
    syntax = 'echo "%s\\t%s" > %s.temp ; bedtools makewindows -g %s.temp -w %s > %s ; rm -rf %s.temp' % (chrom, chrom_length, chrom, chrom, bin_size, file_save_name, chrom)
    process = subprocess.Popen(syntax, stdin=subprocess.PIPE, shell=True,)
    process.communicate()
    if process.returncode != 0:
        sys.stderr.write("Error: Cannot bedtools MAKEWINDOWS.\nCheck if bedtools2 (https://github.com/arq5x/bedtools2) has been installed correctly.\n")
        sys.exit(1)


def main(args, pool):
    chroms = args.species_chroms
    chrom_lengths = args.species_chrom_lengths

    list_of_args = []
    for i, chrom in enumerate(chroms):
        if chrom in chrom_lengths.keys():
            chrom_length = chrom_lengths[chrom]
        else:
            print("Can not find the length of ", chrom)
        list_of_args.append((chrom, chrom_length, args.bin_size))

    create_bed_windows_partial = partial(create_bed_windows)
    pool.starmap(create_bed_windows_partial, list_of_args)
