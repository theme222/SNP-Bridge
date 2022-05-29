import numpy as np
from random import randint, choice, shuffle
import matplotlib.pyplot as plt
from TechnicalTools import log, time_convert, multi_delete, color_gen
from time import time
import json


class SNP:
    """ class to store SNP data type containing allele and global position"""

    def __init__(self, start_pos, local=None):
        self.allele = []
        self.global_start_pos = start_pos
        if local is None:
            self.local_start_pos = start_pos
        else:
            self.local_start_pos = local

    def __repr__(self):
        return f'{self.global_start_pos}'

    '''
    def get_allele(self, *args):
        if not self.allele:  # self.allele == []
            for i in args:
                self.allele.append(i.read[self.start_pos])
        return self.allele
    '''

    def snp_checker(self, snplist):  # check whether a snp is in a lsist
        try:
            index = [i.global_start_pos for i in snplist].index(self.global_start_pos)
            return True, index
        except ValueError:
            return False, None

    @staticmethod
    def key(snp):
        return snp.global_start_pos


class DNA:
    """
    class to store DNA type information
    """

    nucleotides = ['A', 'T', 'G', 'C']
    current_insertid = 1  # keeping track and making sure insert Id aren't the same
    global_snp = []  # all snps in the current DNA strand

    def __init__(self, read=np.array([], dtype=str), snp=None, start_pos=0):
        if snp is None:
            snp = []
        self.read = read  # Numpy array
        self.snp = snp  # Normal Python list
        self.start_pos = start_pos  # Int
        self.insertid = 0  # Int

    def __repr__(self):
        return f'"{self.start_pos}"'
        #  return f'{list(self.read)} \n mutates at {self.snp}'

    def simulate(self, size):
        self.read = np.random.choice(DNA.nucleotides, size + 1)

    '''
    def mutate(self, snp_count):
        for l in range(snp_count):
            index = randint(1, self.size) - 1
            self.read[index] = choice(DNA.nucleotides)
            self.snp = np.append(self.snp, index)
        self.snp.sort()
    '''

    @staticmethod
    def mutate(type_a, type_b, snp_count):
        snp_pos = []
        for i in range(snp_count):
            r = randint(1, type_a.size) - 1
            if r not in snp_pos:
                snp_pos.append(r)

        for l in snp_pos:
            current_type = choice([type_a, type_b])
            # making sure the choice isn't the same as the original
            DNA.nucleotides.remove(current_type.read[l])
            current_type.read[l] = choice(DNA.nucleotides)
            DNA.nucleotides = ['A', 'T', 'G', 'C']

            current_type.snp.append(SNP(l))
            DNA.global_snp.append(SNP(l))
        type_a.snp.sort(key=SNP.key)
        type_b.snp.sort(key=SNP.key)
        DNA.global_snp.sort(key=SNP.key)

    """
    Adenine == 1
    Thymine == 2
    Cytosine == 3
    Guanine == 4
    """

    def to_int(self):

        if self.read.dtype == int:
            return self.read
        elif self.read.dtype == '<U1':

            if len(self.read.shape) == 1:
                _temp = np.array([], dtype=int)  # return variable
                for acid in self.read:
                    if acid == 'A' or acid == 'a':
                        _temp = np.append(_temp, 1)
                    elif acid == 'T' or acid == 't':
                        _temp = np.append(_temp, 2)
                    elif acid == 'C' or acid == 'c':
                        _temp = np.append(_temp, 3)
                    elif acid == 'G' or acid == 'g':
                        _temp = np.append(_temp, 4)
                    else:
                        log('error', 'Unknown character in to int function', acid)
                return _temp
            else:
                _temp2 = np.zeros((0, self.read.shape[1]), dtype=int)
                for read in self.read:
                    _temp = np.array([], dtype=int)  # return variable
                    for acid in read:
                        if acid == 'A' or acid == 'a':
                            _temp = np.append(_temp, 1)
                        elif acid == 'T' or acid == 't':
                            _temp = np.append(_temp, 2)
                        elif acid == 'C' or acid == 'c':
                            _temp = np.append(_temp, 3)
                        elif acid == 'G' or acid == 'g':
                            _temp = np.append(_temp, 4)
                        else:
                            log('error', 'Unknown character in to int function', acid)
                    _temp2 = np.concatenate((_temp2, [_temp]))
                return _temp2

    def to_str(self):

        if self.read.dtype == '<U1':
            return self.read
        elif self.read.dtype == int:

            if len(self.read.shape) == 1:
                _temp = np.array([], dtype=int)  # return variable
                for acid in self.read:
                    if acid == 1:
                        _temp = np.append(_temp, 'A')
                    elif acid == 2:
                        _temp = np.append(_temp, 'T')
                    elif acid == 3:
                        _temp = np.append(_temp, 'C')
                    elif acid == 4:
                        _temp = np.append(_temp, 'G')
                    else:
                        log('error', 'Unknown character in to asc function', acid)
                return _temp
            else:
                _temp2 = np.zeros((0, self.read.shape[1]), dtype=int)
                for read in self.read:
                    _temp = np.array([], dtype=int)  # return variable
                    for acid in read:
                        if acid == 1:
                            _temp = np.append(_temp, 'A')
                        elif acid == 2:
                            _temp = np.append(_temp, 'T')
                        elif acid == 3:
                            _temp = np.append(_temp, 'C')
                        elif acid == 4:
                            _temp = np.append(_temp, 'G')
                        else:
                            log('error', 'Unknown character in to asc function', acid)
                    _temp2 = np.concatenate((_temp2, [_temp]))
                return _temp2

    @property
    def size(self):
        return len(self.read)

    @classmethod
    def derive(cls, obj):
        return cls(read=obj.read.copy())

    @classmethod
    def obtain_read(cls, obj, read_size, insert_size, index=None, insertid=0):  # Makes one insert
        if index is None:
            index = randint(0, obj.size - insert_size)

        # Goes through the global snp list and finds a SNP object to then make a new object with global and local position
        forwardSNP = [SNP(i.global_start_pos, local=i.global_start_pos - index) for i in obj.global_snp
                      if index < i.global_start_pos < index + read_size]

        reverseSNP = [SNP(i.global_start_pos, local=i.global_start_pos - index - insert_size + read_size) for i in
                      obj.global_snp
                      if index + insert_size - read_size < i.global_start_pos < index + insert_size]

        forward = cls(read=obj.read[index:index + read_size], snp=forwardSNP, start_pos=index)
        reverse = cls(read=obj.read[index + insert_size - read_size:index + insert_size], snp=reverseSNP,
                      start_pos=index + insert_size - read_size)
        forward.insertid = insertid
        reverse.insertid = insertid
        return [forward, reverse]

    @staticmethod
    def blend(strand, avg_depth, read_size, insert_size):  # Makes multiple inserts
        return_list = np.array([], dtype=object)
        insert_start_pos = np.random.randint(0, strand.size - insert_size,
                                             round((strand.size * avg_depth) / read_size / 2))
        insert_start_pos.sort()
        for i in insert_start_pos:
            return_list = np.append(return_list, DNA.obtain_read(
                strand, read_size, insert_size, index=i, insertid=DNA.current_insertid))
            DNA.current_insertid += 1
        return return_list

    @staticmethod
    def find_rev_by_insertid(read, shuffled_reads):
        insertid = read.insertid
        for i in shuffled_reads:
            if i.insertid == insertid and i != read:
                return i


def display_simulation(
        template, type_a, type_b, blend_a, blend_b, bridge):  # just a way to plot the results of the simulation
    y_graph = 0  # A standard variable determining what the y value for every graph so it doesn't go on top of each other

    plt.style.use('seaborn-whitegrid')
    fig, (plot1, plot2) = plt.subplots(figsize=(10, 10), nrows=2, ncols=1)
    plot1.set_title(f'Computer DNA read simulation size : {blend_a.size}')
    plot2.set_title(f'Bridge results with {len(bridge)} groups')

    # plotting baseline
    plot1.plot([0, template.size], np.full(2, y_graph), 'k-')
    plot2.plot([0, template.size], np.full(2, y_graph), 'k-')
    y_graph += 1

    ''' ============================================== Making plot 1 ============================================== '''

    # plotting chopped sizes

    _temp2 = y_graph
    _temp3 = 0

    afor = None
    arev = None
    bfor = None
    brev = None

    for i, v in enumerate(blend_a):  # plotting strands from blend_a
        if i % 2 == 0:
            _temp3 = v.start_pos + v.size
            afor, = plot1.plot([v.start_pos, v.start_pos + v.size], np.full(2, y_graph), 'c-')
        else:
            arev, = plot1.plot([v.start_pos, v.start_pos + v.size], np.full(2, y_graph), 'b-')
            plot1.plot([_temp3, v.start_pos], np.full(2, y_graph), 'g--')
            y_graph += 1

    _temp2 += 1

    for i, v in enumerate(blend_b):  # plotting strands from blend_b
        if i % 2 == 0:
            _temp3 = v.start_pos + v.size
            bfor, = plot1.plot([v.start_pos, v.start_pos + v.size], np.full(2, y_graph), 'y-')
        else:
            brev, = plot1.plot([v.start_pos, v.start_pos + v.size], np.full(2, y_graph), 'r-')
            plot1.plot([_temp3, v.start_pos], np.full(2, y_graph), 'g--')
            y_graph += 1

    _temp2 += 5

    # plotting differences
    y_graph += 1
    snp, = plot1.plot([i.global_start_pos for i in DNA.global_snp], np.full(len(DNA.global_snp), y_graph), 'mx',
                      label='Type A SNPS')
    plot1.legend(handles=[afor, arev, bfor, brev, snp],
                 labels=['A Forward Reads', 'A Reverse Reads', 'B Forward Reads', 'B Reverse Reads', 'SNPS'],
                 loc=6)

    ''' ============================================== Making plot 2 ============================================== '''

    for index, group in enumerate(bridge):
        color_ = color_gen('name')
        for read in group:
            plot2.plot([read.start_pos, read.start_pos + read.size], np.full(2, index + 1), color_)

    plt.show()


def bridge_strands(shuffled_reads, snplist):  # function independent from simulator
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

        # find reverse
        _rev1 = []
        for read in group1:
            _rev1.append(DNA.find_rev_by_insertid(read, shuffled_reads))
        group1.extend(_rev1)

        _rev2 = []
        for read in group2:
            _rev2.append(DNA.find_rev_by_insertid(read, shuffled_reads))
        group2.extend(_rev2)

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
    with open("config.json") as json_data:
        configs = json.load(json_data)

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


if __name__ == '__main__':
    main()
