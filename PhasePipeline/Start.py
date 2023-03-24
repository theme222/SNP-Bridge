import os
import json
from helperfile import log
import shutil
from makerfile import *

with open("config.json") as json_data:
    configs = json.load(json_data)
if configs["Input Filename"].split('.')[-1] == 'fastq':
    log("debug","Activating BWA")
    os.system(f"bwa index {configs['Reference Filename']}")
    os.system(f"bwa mem {configs['Reference Filename']} {configs['Input Filename']} > bwa-output.sam")
    log("debug","Running bridger")
    timer = time()
    blended_reads, reference = sam_reader("bwa-output.sam", configs["Reference Filename"])
    DNA.snp_purger(blended_reads)
    bridge_output = bridge_strands(blended_reads, DNA.global_snp)
    print(bridge_output)
    bridge_output = drill_down(bridge_output)
    log('Results', f'Final bridging simulation with a total of {len(bridge_output)} groups.\n', bridge_output)
    log('debug', 'Total time to simulate :', time_convert(time() - timer))
    if input("Show graph? [y/n] : ") == 'y':
        display_simulation(len(reference), bridge_output)
    if input("Write to fasta? [y/n] : ") == 'y':
        fasta_writer(bridge_output, "bwa-output.sam", configs["Reference Filename"], configs['Output Filename'])

    log("debug", "Removing reference index files")
    try:
        os.remove(f"{configs['Reference Filename']}.amb")
        os.remove(f"{configs['Reference Filename']}.ann")
        os.remove(f"{configs['Reference Filename']}.bwt")
        os.remove(f"{configs['Reference Filename']}.pac")
        os.remove(f"{configs['Reference Filename']}.sa")
        os.remove("bwa-output.sam")
    except FileNotFoundError:
        log("error","File not found")