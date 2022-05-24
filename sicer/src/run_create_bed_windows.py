from functools import partial
import subprocess

def create_bed_windows(chrom, chrom_length, bin):
    file_save_name = chrom + '.windows'
    syntax = 'echo -e "%s\\t%s" > %s.temp ; bedtools makewindows -g %s.temp -w %s > %s ; rm -rf %s.temp' % (chrom, chrom_length, chrom, chrom, bin, file_save_name, chrom)
    print (syntax)
    subprocess.Popen(syntax, stdin=subprocess.PIPE, shell=True,)

def separate_bedpe_chroms(file, chrom):
    file_name = os.path.basename(file)
    file_name = file_name.replace('.bed', '')

    match = chrom + "[[:space:]]"
    new_file_name = file_name + '_' + chrom + '-results.bed'
    syntax = 'grep %s %s > %s' % (match, file, new_file_name)
    subprocess.Popen(syntax, stdin=subprocess.PIPE, shell=True,)

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

    separate_chroms_partial = partial(separate_bedpe_chroms, file)
    pool.map(separate_chroms_partial, chroms)

