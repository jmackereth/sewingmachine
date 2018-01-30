import numpy as np
import os

class Linelist(object):
    '''
    Linelist
    ----------------
    a object to contain a parsed linelist, or for making linelists (future)
    ----------------
    INPUT:
    file - the linelist file
    ----------------
    HISTORY:
    2018/29/01 - Written - Mackereth (ARI, LJMU)
    '''
    def __init__(self, file):
        if not os.path.isfile(file):
            raise FileError('Linelist file does not exist!')
        self.labels, self.integration, self.windows = parseLinelist(file)


def parseLinelist(file):
    '''
    function for parsing linelist files (as specified in docs)
    '''
    linelist =np.genfromtxt(file, dtype=None, names=True)
    int_reg = []
    cont_reg = []
    for i in range(0,len(linelist)):
        int_tuple = (linelist['i_b'][i], linelist['i_r'][i])
        cont_tuples = eval(linelist['cont'][i])
        int_reg.append(int_tuple)
        cont_reg.append(cont_tuples)
    labels = linelist['Label']
    return labels, int_reg, cont_reg


class Error(Exception):
    pass

class FileError(Error):
    def __init__(self, message):
        self.message = message
