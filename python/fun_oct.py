import numpy as np


def import_SU_dat(filename, size=(1024, 256, 256), mode='complex64'):

    F = open(filename, 'rb')

    A = F.read()
    B = np.fromstring(A, dtype=mode)
    print(B.shape)
    B = np.reshape(B, size)
    F.close()

    return(B)