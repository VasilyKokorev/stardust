from astropy.io import fits
from astropy.table import Table
import numpy as np
import scipy as sc
from scipy import stats
from scipy.stats import norm
from astropy  import units as u
from astropy import constants as const
import scipy as sp
from scipy.integrate import quad
from scipy import integrate
import math
from scipy.signal import convolve
from scipy import interpolate
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from scipy.stats import norm
from matplotlib.pyplot import *

class Filt:
    def __init__(self, name=None, wave=None, throughput=None):
        self.name = name
        self.wave = wave
        self.throughput = throughput

file_path='filters/filters.txt'
with open(file_path, 'r') as fp:
    lines = fp.readlines()

filters = []
wave = []
for line in lines:
    if 'lambda_c' in line:
        if len(wave) > 0:
            new_filter = Filt(name=header,wave=np.cast[float](wave),throughput=np.cast[float](through))
            filters.append(new_filter)
        header = ' '.join(line.split()[1:])
        wave = []
        through = []
    else:
        lspl = np.cast[float](line.split())
        wave.append(lspl[1])
        through.append(lspl[2])

new_filter = Filt(name=header,wave=np.cast[float](wave), throughput=np.cast[float](through))

filters.append(new_filter)

print('Filters Imported')
