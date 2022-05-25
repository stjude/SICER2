# Developed by Zang Lab at University of Virginia - 2018

#Author: Jin Yong Yoo

import os
import shutil
import sys
import tempfile
import multiprocessing as mp

curr_path = os.getcwd()

# From SICER Package
from sicer.src import remove_redundant_reads
from sicer.src import run_make_graph_file_by_chrom
from sicer.src import create_bed_windows
from sicer.src import separate_bedpe_chroms
from sicer.src import process_and_clean_bedpe
from sicer.src import import_graph_file_by_chrom
from sicer.src import find_islands_in_pr
from sicer.src import associate_tags_with_chip_and_control_w_fc_q
from sicer.src import filter_islands_by_significance
from sicer.src import make_normalized_wig
from sicer.src import filter_raw_tags_by_islands

''' args: ArgumentParser object formed form command line parameters
    df_run: If df_run is true, then this instance of SICER is called by SICER-DF module.
            Default value is False.'''


def main(args, df_run=False):
    # Checks if there is a control library
    control_lib_exists = True
    if args.control_file is None:
        control_lib_exists = False

    # Creates temporary directory to contain all intermediate files.
    try:
        temp_dir = tempfile.mkdtemp()
        # Change current working directory to temp_dir
        os.chdir(temp_dir)
        print(temp_dir)

    except:
        sys.exit(
            "Temporary directory required for SICER cannot be created. Check if directories can be created in %s." % curr_path)

    try:
        # Step 0: create Pool object for parallel-Processing
        num_chroms = len(args.species_chroms)
        pool = mp.Pool(processes=min(args.cpu, num_chroms))

        if args.paired_end == True:

            # Step 1-PE: creating bed windows
            print("Creating bed windows... \n")
            create_bed_windows.main(args, pool) #make windows

            # Step 2-PE: Converting the bedpe to graph windows
            print("Create graph bin based on pre-defined bin size %s (bp)... \n" % args.bin_size)
            treatment_file_name = os.path.basename(args.treatment_file)
            separate_bedpe_chroms.main(args, args.treatment_file, pool)
            print('\n')

            print("Partition the genome into bins and create graph files... \n")
            total_tag_in_windows = process_and_clean_bedpe.main(args, args.treatment_file, pool) #bedpe to graph
            args.treatment_file = treatment_file_name
            total_treatment_read_count = total_tag_in_windows
            print('\n')

            # Step 3-PE: Using the control graph file
            if control_lib_exists:
                control_file_name = os.path.basename(args.control_file)
                print("Use the", control_file_name, "control graph file... \n")
                separate_bedpe_chroms.main(args, args.treatment_file, pool)

                total_control_read_count = process_and_clean_bedpe.main(args, args.control_file, pool)
                args.control_file = control_file_name
                print('\n')

        else:
            # Step 1-SE: Remove redundancy reads in input file according to input threshold
            # Output is the total number of reads retained. Represents size of library.
            treatment_file_name = os.path.basename(args.treatment_file)
            print("Preprocess the", treatment_file_name, "file to remove redundancy with threshold of",
                  args.redundancy_threshold, "\n")
            total_treatment_read_count = remove_redundant_reads.main(args, args.treatment_file, pool)
            args.treatment_file = treatment_file_name
            print('\n')

            # Step 2-SE: Remove redundancy reads in control library according to input threshold
            if control_lib_exists:
                control_file_name = os.path.basename(args.control_file)
                print("Preprocess the", control_file_name, "file to remove redundancy with threshold of",
                      args.redundancy_threshold, "\n")
                total_control_read_count = remove_redundant_reads.main(args, args.control_file, pool)
                args.control_file = control_file_name
                print('\n')

            # Step 3-SE: Partition the genome in windows and generate graph files for each chromsome
            print("Partition the genome in windows and generate summary files... \n")
            total_tag_in_windows = run_make_graph_file_by_chrom.main(args, pool)
            print('\n')

        # Step 4: Normalize and generate WIG file
        print("Normalizing graphs by total island filitered reads per million and generating summary WIG file...\n")
        output_WIG_name = (treatment_file_name.replace('.bed', '') + "-W" + str(args.window_size) + "-normalized.wig")
        make_normalized_wig.main(args, output_WIG_name, pool)

        # Step 5: Find candidate islands exhibiting clustering
        print("Finding candidate islands exhibiting clustering...\n")
        find_islands_in_pr.main(args, total_tag_in_windows, pool)
        print("\n")

        # Running SICER with a control library
        if control_lib_exists:
            # Step 6
            print("Calculating significance of candidate islands using the control library... \n")
            associate_tags_with_chip_and_control_w_fc_q.main(args, total_treatment_read_count, total_control_read_count, pool)

            # Step 7: Filter out any significant islands whose pvalue is greater than the false discovery rate
            print("Identify significant islands using FDR criterion\n")
            significant_read_count = filter_islands_by_significance.main(args, 7, pool)  # 7 represents the ith column we want to filtered by
            print("Out of the ", total_treatment_read_count, " reads in ", treatment_file_name, ", ", significant_read_count, " reads are in significant islands\n")

        # Optional Outputs
        if args.significant_reads:
            # Step 8: Filter treatment reads by the significant islands found from step 8
            print("Filtering reads with identified significant islands...\n")
            filter_raw_tags_by_islands.main(args, pool)

            # Step 9: Produce graph file based on the filtered reads from step 9
            print("Making summary graph with filtered reads...\n")
            run_make_graph_file_by_chrom.main(args, pool, True)
            # Step 10: Produce Normalized WIG file
            print("\nNormalizing graphs by total island filitered reads per million and generating summary WIG file...\n")
            output_WIG_name = (treatment_file_name.replace('.bed', '') + "-W" + str(args.window_size) + "-G" + str(args.gap_size) + "-FDR" + str(args.false_discovery_rate) + "-islandfiltered-normalized.wig")
            make_normalized_wig.main(args, output_WIG_name, pool)

        pool.close()
        pool.join()

        # Final Step
        if df_run == True:
            return temp_dir, total_treatment_read_count
        else:
            print("End of SICER")
    finally:
        if df_run == False:
            print("Removing temporary directory and all files in it.")
            #shutil.rmtree(temp_dir)

