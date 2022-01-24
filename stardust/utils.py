import numpy as np


def find_nearest(array, value):
    '''
    Fing the element of the array which is closest to value
    '''
    n = [abs(i-value) for i in array]
    idx = n.index(min(n))
    return idx

def is_number(s):
    '''
    Check if string has numbers
    '''
    try:
        float(s)
        return True

    except ValueError:
        return False

def to_sigma(arr):
        '''
        Compute the uncertainty from confidence intervals
        '''
        try:
            perc = np.percentile(arr,q=[16,50,84])
        except:
            return -99.
        return np.mean([perc[2]-perc[1],perc[1]-perc[0]])
