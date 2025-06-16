import numpy as np
from astropy.time import Time
import math

def which_detectors(detectors_string, detectors=["H1","L1","V1"]):
    detectors_list=[]
    for i in detectors:
        if i in detectors_string:
            detectors_list.append(i)
    return detectors_list

def gps2jd(timeGPS):
    t=Time(timeGPS, format="gps")
    t.format="jd"
    #print(t.value)
    return t.value
