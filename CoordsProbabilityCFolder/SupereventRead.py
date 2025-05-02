from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from hpmoc import PartialUniqSkymap
def superevent_read(event_array, file_location="./O4a_fits/"):
    superevent_fits=file_location+event_array[0]+"_bayestar_fits_gz_"+str(event_array[1])+".fits"
    return PartialUniqSkymap.read(superevent_fits, strategy="ligo")
