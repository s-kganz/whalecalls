import numpy as np 
import warnings
from functools import reduce
import operator

def find_longest_sequence(arr):
    '''Finds longest consecutive sequence in array.'''
    start = 0
    seq = [0, 0]
    for i in range(0, len(arr) - 1):
        if arr[i + 1] - arr[i] != 1 and i + 1 - start > len(seq):
            seq = arr[start:i+1]
            start = i + 1
    return seq[0], seq[-1]

def sub_avg(s):
    '''Subtract mean of list from each element of the list.'''
    new = s.copy()
    means = np.mean(new, axis=1).reshape(len(new), 1)
    return new - means

def group_consecutives(vals, step=1):
    '''Returns 2d list containing all consecutive subsets of array based off of given step size.'''
    run = []
    result = [run]
    expect = None
    for v in vals:
        if v == expect or expect is None:
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step
    return result


