import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const



"""
@verison: 1

@author: David Wilson

@date 20190809

Turns XMM data in to HLSP format ecsv and fits files, ready to be added to an SED
"""

