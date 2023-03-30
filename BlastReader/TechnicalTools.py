"""
This is the Technical tools python library created for general purpose stuff
just so you know I have no idea what half of this means so don't ask
"""

from inspect import currentframe, getframeinfo
from random import randint, choice
import os


def log(logtype, message, var=None):
    """
    like the logging import but less ugly
    :param var:
    :param message:
    :param logtype:
    """
    linenumber = currentframe().f_back.f_lineno
    if os.name == 'nt':
        filename = getframeinfo(currentframe().f_back).filename.split('\\')[-1]
    else:
        filename = getframeinfo(currentframe().f_back).filename.split('/')[-1]
    if logtype == 'info':
        print_message = f'\033[1;34m[{logtype.upper()}] in {filename} at line {linenumber} - {message} |{var}|\033[0m'
    elif logtype == 'debug':
        print_message = f'\033[1;35m[{logtype.upper()}] in {filename} at line {linenumber} - {message} |{var}|\033[0m'
    elif logtype == 'error':
        print_message = f'\033[0;31m[{logtype.upper()}] in {filename} at line {linenumber} - {message} |{var}|\033[0m'
    elif logtype == 'critical':
        print_message = f'\033[1;31m[{logtype.upper()}] in {filename} at line {linenumber} - {message} |{var}|\033[0m'
    else:
        print_message = f'\033[32m[{logtype.upper()}] in {filename} at line {linenumber} - {message} |{var}|\033[0m'
    if var is None:
        print(print_message[:-10])
    else:
        print(print_message)


def time_convert(seconds):
    """
    turns absurd amounts of seconds into understandable text
    :param seconds:
    :return:
    """

    year = 31557600
    month = 2629800
    week = 604800
    day = 86400
    hour = 3600
    minute = 60
    second = 1
    millisecond = 0.001
    microsecond = 0.000001
    nanosecond = 0.000000001

    return_text = ''

    amount = seconds // year
    if amount > 0:
        return_text = return_text + f' {int(amount)} Year'
        seconds -= year * amount

    amount = seconds // month
    if amount > 0:
        return_text = return_text + f' {int(amount)} Month'
        seconds -= month * amount

    amount = seconds // week
    if amount > 0:
        return_text = return_text + f' {int(amount)} Week'
        seconds -= week * amount

    amount = seconds // day
    if amount > 0:
        return_text = return_text + f' {int(amount)} Day'
        seconds -= day * amount

    amount = seconds // hour
    if amount > 0:
        return_text = return_text + f' {int(amount)} Hour'
        seconds -= hour * amount

    amount = seconds // minute
    if amount > 0:
        return_text = return_text + f' {int(amount)} Minute'
        seconds -= minute * amount

    amount = seconds // second
    if amount > 0:
        return_text = return_text + f' {int(amount)} Second'
        seconds -= second * amount

    amount = int(seconds / millisecond)
    if amount > 0:
        return_text = return_text + f' {amount} MilliSecond'
        seconds -= millisecond * amount

    amount = int(seconds / microsecond)
    if amount > 0:
        return_text = return_text + f' {amount} MicroSecond'
        seconds -= microsecond * amount

    amount = int(seconds / nanosecond)
    if amount > 0:
        return_text = return_text + f' {amount} NanoSecond'
        seconds -= nanosecond * amount

    return return_text + ' '


def color_gen(style='rgb'):
    """ Makes a random color depending on style """
    color_names = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'cyan', 'black',
                   'magenta', 'yellow', 'aqua', 'purple', 'pink']
    if style == 'rgb' or style == 'RGB':
        return [randint(0, 255), randint(0, 255), randint(0, 255)]
    elif style == 'hex':
        return '#' + format(randint(0, 16777215), 'x')
    elif style == 'name' or style == 'color':
        return choice(color_names)


class TextFile:
    """Simple class to format text files"""

    def __init__(self):
        pass

    def __repr__(self):
        return 'Why would you do this?'

    @staticmethod
    def delnewline(filename, char_amount):
        """Deletes a set number of characters after a new line"""
        finishedtext = ''
        for line in open(filename, 'r'):
            finishedtext += line[char_amount:]
        filename_new = filename.split('.')[0] + '_new.txt'
        try:
            open(filename_new, 'x').write(finishedtext)
        except FileExistsError:
            log('error', 'Created file already exists', filename_new)

    @staticmethod
    def replaceword(filename, word, replacement=''):
        """finds the word and replaces with another thing"""
        finishedtext = ''
        for line in open(filename, 'r'):
            finishedtext += line.replace(word, replacement)
        filename_new = filename.split('.')[0] + '_new.txt'
        try:
            open(filename_new, 'x').write(finishedtext)
        except FileExistsError:
            log('error', 'Created file already exists', filename_new)


def multi_delete(list_, args):
    args = set(args)
    indexes = sorted(args, reverse=True)
    for index in indexes:
        del list_[index]
    return list_


def main():
    log('critical', 'not critical just epic')
    log('debug', 'Time test thing lol', time_convert(604860.1238349))


if __name__ == '__main__':
    main()
