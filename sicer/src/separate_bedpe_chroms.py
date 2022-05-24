from functools import partial
import subprocess

def separate_bedpe_chroms(file, chrom):
    file_name = os.path.basename(file)
    file_name = file_name.replace('.bed', '')

    match = chrom + "[[:space:]]"
    new_file_name = file_name + '_' + chrom + '-results.bed'
    syntax = 'grep %s %s > %s' % (match, file, new_file_name)
    subprocess.Popen(syntax, stdin=subprocess.PIPE, shell=True,)

def main(args, file, pool):
    chroms = args.species_chroms

    separate_chroms_partial = partial(separate_bedpe_chroms, file)
    pool.map(separate_chroms_partial, chroms)
