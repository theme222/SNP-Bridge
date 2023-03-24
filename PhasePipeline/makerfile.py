import matplotlib.pyplot as plt
from time import time
from helperfile import *
import json


def display_simulation(size, bridge):  # just a way to plot the results of the simulation
    y_graph = 0  # A variable determining what the y value for every graph, so it doesn't go on top of each other

    plt.style.use('seaborn-whitegrid')
    fig, plot2 = plt.subplots()
    plot2.set_title(f'Final phase created : {len(bridge)}')

    # plotting baseline

    plot2.plot([0, size], np.full(2, y_graph), 'k-')
    y_graph += 1

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
    for index, bridge in enumerate(return_list):
        bridge.sort(key=DNA.key)
        return_list[index] = bridge
    return return_list


def fasta_writer(groups, samfile, ref,fn):
    list_txt = []
    for group in groups:
        text = group[0].read_with_ins
        starter = group[0].start_pos
        for read in group[1:]:
            text = text[:read.start_pos - starter] + read.read_with_ins
        list_txt.append([text,starter])
    with open(fn, 'w') as fasta_file:
        for index, txt in enumerate(list_txt):
            fasta_file.write(f">SamSplicer.py-{index} Phase output from bridging with {samfile} using reference {ref}"
                             + f" length {len(txt[0])} start-pos {txt[1]} \n")
            fasta_file.write(txt[0] + '\n')


def main():
    with open("config.json") as json_data:
        configs = json.load(json_data)
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
        fasta_writer(bridge_output, "bwa-output.sam", configs["Reference Filename"],configs['Output Filename'])


if __name__ == '__main__':
    main()
