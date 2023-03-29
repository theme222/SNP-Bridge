import os
import json
from helperfile import log
import shutil
from makerfile import *


# ------------- AI Code Incoming -------------
def get_file_list(src_dir):
    # Get a list of all files in the subdirectory
    file_list = os.listdir(src_dir)
    return file_list


def copy_file(src_dir, dest_dir, file_name, file_index):
    # Get the list of files from the separate function
    # Copy the next file to the main directory and rename it
    full_file_name = os.path.join(src_dir, file_name)
    if os.path.isfile(full_file_name):
        new_file_name = f'reference{file_index + 1}.fasta'
        shutil.copy(full_file_name, os.path.join(dest_dir, new_file_name))
        # Returns the new file name
        return new_file_name
    else:
        raise FileNotFoundError


# ------------- AI Code Exiting -------------


def main():
    with open("config.json") as json_data:
        configs = json.load(json_data)
    src_dir = 'References'
    dest_dir = '.'
    reference_list = get_file_list(src_dir)
    best_output = [i for i in range(69)]
    best_file = ''

    for index, file_name in enumerate(reference_list):
        log('info', "Using reference", file_name)

        current_reference_filename = copy_file(src_dir, dest_dir, file_name, index)

        if configs["Input Filename"].split('.')[-1] == 'fastq':
            log("debug", "Activating BWA")
            os.system(f"bwa index {current_reference_filename}")
            os.system(f"bwa mem {current_reference_filename} {configs['Input Filename']} > bwa-output.sam")
            log("debug", "Running bridger")
            timer = time()
            blended_reads, reference = sam_reader("bwa-output.sam", current_reference_filename)
            DNA.snp_purger(blended_reads)
            bridge_output = bridge_strands(blended_reads, DNA.global_snp)
            print(bridge_output)
            bridge_output = drill_down(bridge_output)
            if configs["Output Type"] == 'Best' and len(bridge_output) < len(best_output):
                best_output = bridge_output
                best_file = file_name
            log('Results',
                f'Final bridging simulation with a total of {len(bridge_output)} groups. Using reference {file_name}\n',
                bridge_output)
            log('debug', 'Total time to simulate (with BWA) :', time_convert(time() - timer))
            try:
                if configs["Output Type"] == 'Select':
                    if input("Write to FastA? [y/n] : ") == 'y':
                        fasta_writer(bridge_output, 'bwa-output.sam', file_name)
                elif configs["Output Type"] == "All":
                    fasta_writer(bridge_output, 'bwa-output.sam', file_name)
            except KeyboardInterrupt:
                log('debug', "KeybardInterrupt")

            log("debug", "Removing reference index files")

            try:
                os.remove(f"{current_reference_filename}.amb")
                os.remove(f"{current_reference_filename}.ann")
                os.remove(f"{current_reference_filename}.bwt")
                os.remove(f"{current_reference_filename}.pac")
                os.remove(f"{current_reference_filename}.sa")
                os.remove(f"{current_reference_filename}")
                os.remove("bwa-output.sam")
            except FileNotFoundError:
                log("error", "File not found")

        if best_output != [i for i in range(69)]:
            fasta_writer(best_output, 'bwa-output.sam', best_file)


if __name__ == '__main__':
    main()
