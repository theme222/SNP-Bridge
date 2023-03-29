import numpy as np
from random import randint, choice, shuffle
import matplotlib.pyplot as plt
from TechnicalTools import log, time_convert, multi_delete, color_gen
from time import time
from SamReader import *
import json
from Blaster import blast
from TechnicalTools import log
import threading


def display_simulation(size, bridge):  # just a way to plot the results of the simulation
    y_graph = 1  # A variable determining what the y value for every graph, so it doesn't go on top of each other

    plt.style.use('seaborn-whitegrid')
    fig, plot1 = plt.subplots()
    plot1.set_title(f'Final phase created : {len(bridge)}')

    # plotting baseline

    ref, = plot1.plot([0, size], np.full(2, 0), 'k-')

    _temp = []
    for group in (bridge):
        color_ = color_gen('name')
        _temp2 = None
        for read in group:
            _temp2, = plot1.plot([read.start_pos, read.start_pos + read.size], np.full(2, y_graph), color_)
        _temp.append(_temp2)
        y_graph += 1

    snp, = plot1.plot([i.global_start_pos for i in DNA.global_snp], np.full(len(DNA.global_snp), y_graph), 'mx',
                      label='Type A SNPS')
    _temp3 = [ref, snp]
    _temp3.extend(_temp)
    _temp4 = ['Reference Seq', 'SNPs']
    _temp4.extend([f"Phase {i}" for i in range(len(_temp3))])
    plot1.legend(handles=_temp3,
                 labels=_temp4,
                 loc=6)
    plt.show()


def bridge_strands(shuffled_reads, snplist):  # function independent of simulator
    return_list = []
    bridge_set1 = []
    bridge_set2 = []
    stable1 = False
    stable2 = False
    for snp in snplist:

        reads_within_snp = []

        # find the reads that contain the snp
        for read in shuffled_reads:
            contains_snp, index = snp.snp_checker(read.snp)
            if contains_snp:
                reads_within_snp.append([read, index])

        # separate the reads into two groups

        pivot = None
        group1 = []
        group2 = []
        garbage = []
        for readset in reads_within_snp:
            if pivot is None and readset[0].read[readset[0].snp[readset[1]].local_start_pos] in snp.values:
                pivot = readset[0].read[readset[0].snp[readset[1]].local_start_pos]
                group1.append(readset[0])
            else:
                if readset[0].read[readset[0].snp[readset[1]].local_start_pos] == pivot:
                    group1.append(readset[0])
                elif readset[0].read[readset[0].snp[readset[1]].local_start_pos] in snp.values:
                    group2.append(readset[0])
                else:
                    garbage.append(readset[0])

        # checking if the set is empty
        if not bridge_set1 and group1:  # bridge set == []
            bridge_set1.append(group1)
        if not bridge_set2 and group2:
            bridge_set2.append(group2)

        if bridge_set1 and bridge_set2:  # precaution for Index Error
            switched = False
            for i in group1:  # Checking if we need to switch the two groups
                if i in bridge_set2[-1]:
                    group1, group2 = group2, group1
                    switched = True
                    break
            if not switched:  # performance optimization
                for i in group2:  # Checking if we need to switch the two groups
                    if i in bridge_set1[-1]:
                        group1, group2 = group2, group1
                        break

        # checking connection between snps
        if bridge_set1:  # precaution for Index Error
            for i in group1:
                #  log('debug', f'checking if {i} is in {bridge_set1[-1]} g1')
                if i in bridge_set1[-1]:
                    log('info', 'SNP Connected g1', f'{snp}')
                    stable1 = True
                    bridge_set1.append(group1)
                    break

        if bridge_set2:  # precaution for Index Error
            for i in group2:
                #  log('debug', f'checking if {i} is in {bridge_set2[-1]} g2')
                if i in bridge_set2[-1]:
                    log('info', 'SNP Connected g2', f'{snp}')
                    stable2 = True
                    bridge_set2.append(group2)
                    break

        if not stable1 and bridge_set1:
            log('error', 'Unsuccessful SNP connection g1', f'{snp}')
            return_list.append(bridge_set1)
            if group1:
                bridge_set1 = [group1]
            else:
                bridge_set1 = []
        stable1 = False

        if not stable2 and bridge_set2:
            log('error', 'Unsuccessful SNP connection g2', f'{snp}')
            return_list.append(bridge_set2)
            if group2:
                bridge_set2 = [group2]
            else:
                bridge_set2 = []
        stable2 = False

    if bridge_set1:
        return_list.append(bridge_set1)
    if bridge_set2:
        return_list.append(bridge_set2)
    return return_list


def drill_down(chonky_list):
    return_list = []
    for group in chonky_list:
        if group:
            tiny_return_list = []
            for tiny_group in group:
                for read in tiny_group:
                    if read not in tiny_return_list:
                        tiny_return_list.append(read)
            return_list.append(tiny_return_list)

    for index, bridge in enumerate(return_list):
        bridge.sort(key=DNA.key)
        return_list[index] = bridge
    return return_list


def fasta_maker(groups, samfile, ref,ifblast=True):
    blast_thread = None
    list_txt = []
    for group in groups:
        # Gather the insertions in the group
        insertions = np.array([])
        for read in group:  # AI code to remove the duplicate insertions
            insertions = np.append(insertions, read.insertions)
            _, indices = np.unique([insertion.global_start_pos for insertion in insertions], return_index=True)
            insertions = insertions[np.sort(indices)[::-1]]  # sort from most to least

        text = group[0].read
        starter = group[0].start_pos
        for read in group[1:]:  # create read without insertions
            text = np.append(text[:read.start_pos - starter], read.read)
        # insert the insertions 😂
        for ins in insertions:
            text = np.insert(text, ins.local_start_pos, ins.values)

        list_txt.append(["".join(text), starter])  # sequence, index
    if ifblast:
        log('info', 'Beginning Thread')
        blast_thread = threading.Thread(
            target=blast,
            args=([list(zip([seq[0] for seq in list_txt],[f'phase{i}'for i in range(len(list_txt))]))])
            )
        blast_thread.start()
    with open("final_output.fasta", 'w') as fasta_file:
        for index, txt in enumerate(list_txt):
            fasta_file.write(f">SamSplicer.py-{index} Phase output from bridging with {samfile} using reference {ref}"
                             + f" length {len(txt[0])} start-pos {txt[1]} \n")
            fasta_file.write(txt[0] + '\n')
    return blast_thread


def main():
    with open("config.json") as json_data:
        configs = json.load(json_data)
    timer = time()
    blended_reads, reference = sam_reader(configs["Input Sam File"], configs["Reference File"])
    DNA.snp_purger(blended_reads)
    bridge_output = bridge_strands(blended_reads, DNA.global_snp)
    print(bridge_output)
    bridge_output = drill_down(bridge_output)
    log('Results', f'Final bridging simulation with a total of {len(bridge_output)} groups.\n', bridge_output)
    log('debug', 'Total time to simulate :', time_convert(time() - timer))
    thread = None
    if input("Write to fasta? [y/n] : ") == 'y':
        thread = fasta_maker(bridge_output, configs["Input Sam File"], configs["Reference File"])
    if input("Show graph? [y/n] : ") == 'y':
        display_simulation(len(reference), bridge_output)
    if thread is not None:
        thread.join()





if __name__ == '__main__':
    main()
