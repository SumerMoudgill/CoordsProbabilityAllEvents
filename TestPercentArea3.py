import numpy as np
import math
from CoordsProbabilityCFolder.PercentArea import *
from CoordsProbabilityCFolder.O4aSuperevents import *
from CoordsProbabilityCFolder.DetectorMinima import *

for i in range(len(O4a_May_events_1)):
    current_superevent_array=O4a_May_events_1[i]
    print("Event name: ", current_superevent_array[0])
    print("Event time: ", current_superevent_array[2])
    file_location="Downloads/CoordsProbabilityAllEvents/O4a_fits/"
    superevent_fits=file_location+current_superevent_array[0]+"_bayestar_fits_gz_"+str(current_superevent_array[1])+".fits"
    print(superevent_fits)
    print("50% area: ", get_area(file_location, current_superevent_array, 50))
    print("90% area: ", get_area(file_location, current_superevent_array, 90))
    print("50% area distance from coordinate zero: ", get_area_distance(file_location, current_superevent_array, [0,0], 50, order=8))
    print("90% area distance from coordinate zero: ", get_area_distance(file_location, current_superevent_array, [0,0], 90, order=8))
    cf_list=get_all_null_coords(current_superevent_array[2], which_detectors(current_superevent_array[3]))
    print("cf_list:", cf_list)
    res_sep_min_50=np.zeros((3, 4))
    res_sep_min_90=np.zeros((3, 4))
    for i in range(3):
        for j in range(4):
            if len(cf_list[i]) !=0:
                if len(cf_list[i][j]) != 0:
                    current_point=cf_list[i][j]
                    res_sep_min_50[i][j]=get_area_distance(file_location, current_superevent_array, current_point, 50, order=8)
                    res_sep_min_90[i][j]=get_area_distance(file_location, current_superevent_array, current_point, 90, order=8)
    print("angular separations from 90% area: ", res_sep_min_90.tolist())
    print("angular separations from 50% area: ", res_sep_min_50.tolist())
    print("Event name: ", current_superevent_array[0])
    print("RES.OUT")
