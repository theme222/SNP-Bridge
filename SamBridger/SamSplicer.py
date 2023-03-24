import numpy as np
from random import randint, choice, shuffle
import matplotlib.pyplot as plt
from TechnicalTools import log, time_convert, multi_delete, color_gen
from time import time
from SamReader import *
import json




def display_simulation(
        size,  bridge):  # just a way to plot the results of the simulation
    y_graph = 0  # A variable determining what the y value for every graph, so it doesn't go on top of each other

    plt.style.use('seaborn-whitegrid')
    fig, (plot1, plot2) = plt.subplots(figsize=(10, 10), nrows=2, ncols=1)
    plot2.set_title(f'Computer DNA read simulation size : {2}')

    # plotting baseline

    plot2.plot([0, size], np.full(2, y_graph), 'k-')
    y_graph += 1

    ''' ============================================== Making plot 2 ============================================== '''

    for index, group in enumerate(bridge):
        color_ = color_gen('name')
        for read in group:
            plot2.plot([read.start_pos, read.start_pos + read.size], np.full(2, index + 1), color_)

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
        for readset in reads_within_snp:
            if pivot is None:
                pivot = readset[0].read[readset[0].snp[readset[1]].local_start_pos]
                group1.append(readset[0])
            else:
                if readset[0].read[readset[0].snp[readset[1]].local_start_pos] == pivot:
                    group1.append(readset[0])
                else:
                    group2.append(readset[0])
        '''
        # find reverse
        _rev1 = []
        for read in group1:
            _rev1.append(DNA.find_rev_by_insertid(read, shuffled_reads))
        group1.extend(_rev1)

        _rev2 = []
        for read in group2:
            _rev2.append(DNA.find_rev_by_insertid(read, shuffled_reads))
        group2.extend(_rev2)
        '''
        #  print('\ngroup 1:', group1)
        #  print('group 2:', group2)

        # checking if the set is empty
        if not bridge_set1 and group1:  # bridge set == []
            bridge_set1.append(group1)
        elif not bridge_set2 and group2:
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

    print(return_list)
    indexes = []
    len_list = [len(i) for i in return_list]
    is_break = False
    for index, group in enumerate(return_list[:-1]):  # checks everything except last
        for read in group[:len_list[index]]:
            _temp = []
            for i, e in enumerate(return_list[index + 1:]):  # checks everything after the thing its checking
                if read in e:
                    indexes.append(index)
                    for r in group:
                        if r not in e:
                            _temp.append(r)
                    e.extend(_temp)
                    is_break = True
                    break
            if is_break:
                is_break = False
                break

    multi_delete(return_list, indexes)
    return return_list


def main():
    timer = time()
    blended_reads = sam_reader("SamFile.sam","ReferenceFile.fasta")
    DNA.snp_purger(blended_reads)
    bridge_output = bridge_strands(blended_reads, DNA.global_snp)
    print(bridge_output)
    bridge_output = drill_down(bridge_output)
    log('Results', f'Final bridging simulation with a total of {len(bridge_output)} groups.\n', bridge_output)
    log('debug', 'Total time to simulate :', time_convert(time() - timer))
    display_simulation(10340,bridge_output)

    '''
    timer = time()
    with open("config.json") as json_data:
        configs = json.load(json_data)
    sequences = fastA_parser("FastAFile.fasta", 2)

    # returns the changed seq 0 a middle line which depicts the changes and seq 1
    abc = sequenceAligner(sequences[0].seq, sequences[1].seq)
    type_a, type_b = DNA(read=np.array(list(abc[0]))), DNA(np.array(list(abc[2])))  # Changes from text to np array

    for i in range(len(abc[0])):  # Find SNPs in FastA (Indels are replaced with X)
        seq0 = abc[0]
        seq1 = abc[2]

        if seq0[i] != seq1[i]:
            DNA.global_snp.append(SNP(i))
    print(DNA.global_snp)

    #  Randomly Assign smaller reads from a big read (type a , type b)
    blend_a = DNA.blend(type_a, configs['Average Depth'], configs['Read Size'], configs['Insert Size'])
    blend_b = DNA.blend(type_b, configs['Average Depth'], configs['Read Size'], configs['Insert Size'])
    print('blend_a', blend_a)
    print('blend b', blend_b)
    bridge_parts = np.append(blend_a, blend_b)
    np.random.shuffle(bridge_parts)
    bridge_output = bridge_strands(bridge_parts, DNA.global_snp)
    print(bridge_output)
    bridge_output = drill_down(bridge_output)
    log('Results', f'Final bridging simulation with a total of {len(bridge_output)} groups.\n', bridge_output)
    log('info', 'Simulating with these configurations\n', configs)
    log('debug', 'Total time to simulate :', time_convert(time() - timer))
    display_simulation(type_a.size, type_a, type_b, blend_a, blend_b, bridge_output)
'''

'''
    DNA.current_insertid = 1
    DNA.global_snp = []
    template = DNA()
    template.simulate(configs['Template Size'])
    type_a = DNA.derive(template)
    type_b = DNA.derive(template)
    DNA.mutate(type_a, type_b, configs['SNP Count'])
    blend_a = DNA.blend(type_a, configs['Average Depth'], configs['Read Size'], configs['Insert Size'])
    blend_b = DNA.blend(type_b, configs['Average Depth'], configs['Read Size'], configs['Insert Size'])
    print('blend_a', blend_a)
    print('blend b', blend_b)
    bridge_parts = np.append(blend_a, blend_b)
    np.random.shuffle(bridge_parts)
    bridge_output = bridge_strands(bridge_parts, DNA.global_snp)
    print(bridge_output)
    bridge_output = drill_down(bridge_output)
    log('Results', f'Final bridging simulation with a total of {len(bridge_output)} groups.\n', bridge_output)
    log('info', 'Simulating with these configurations\n', configs)
    log('debug', 'Total time to simulate :', time_convert(time() - timer))
    display_simulation(template, type_a, type_b, blend_a, blend_b, bridge_output)
'''

if __name__ == '__main__':
    main()
