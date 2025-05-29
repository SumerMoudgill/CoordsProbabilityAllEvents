import numpy as np
import math
from matplotlib import pyplot as plt
from astropy.table import Table
from hpmoc import PartialUniqSkymap
from hpmoc.plot import get_wcs, plot, gridplot
from hpmoc.points import PointsTuple, Rgba
from CoordsProbabilityCFolder.DetectorMinima import *

def CoordsTuple (coords_array, dispersion=3.0):
    '''
    Converts an array of four [latitude, longitude] coordinates into an astropy PointsTuple.
    '''
    len_c = len(coords_array)
    pts1=PointsTuple(
        [
            (coords_array[0][1], coords_array[0][0], dispersion),
            (coords_array[1][1], coords_array[1][0], dispersion),
            (coords_array[2][1], coords_array[2][0], dispersion),
            (coords_array[3][1], coords_array[3][0], dispersion),
        ]
    )
    return pts1

def CoordsScatter (coords_array, skymap_m, filename, dispersion=3.0):
    '''
    Places a PointsTuple of an array of four [latitude, longitude] coordinates on a skymap.
    '''
    len_c = len(coords_array)
    skymap_m.plot(CoordsTuple(coords_array, dispersion))
    plt.savefig(filename+".png")

