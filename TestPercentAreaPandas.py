import numpy as np
import math
from CoordsProbabilityCFolder.PercentArea import *
from CoordsProbabilityCFolder.O4aSuperevents import *
from CoordsProbabilityCFolder.DetectorMinima import *
from CoordsProbabilityCFolder.ConversionF import *
import pandas as pd

iterator0 = 0
nameslist = []
area_90_list = []
area_50_list = []
sep_90_list = []
sep_50_list = []
print("len(O4a_May_events_1):", len(O4a_May_events_1))
for i in range(len(O4a_May_events_1)):
    print("iterator:", iterator0)
    current_superevent_array=O4a_May_events_1[i]
    print("Event name: ", current_superevent_array[0])
    nameslist.append(current_superevent_array[0])
    print("Event time: ", current_superevent_array[2])
    file_location="Downloads/CoordsProbabilityAllEvents/O4a_fits/"
    superevent_fits=file_location+current_superevent_array[0]+"_bayestar_fits_gz_"+str(current_superevent_array[1])+".fits"
    #print(superevent_fits)
    #print("90% area: ", get_area(file_location, current_superevent_array, 90))
    area_90_list.append(get_area(file_location, current_superevent_array, 90))
    area_50_list.append(get_area(file_location, current_superevent_array, 50))
    #print("50% area: ", get_area(file_location, current_superevent_array, 50))
    #print("50% area distance from coordinate zero: ", get_area_distance(file_location, current_superevent_array, [0,0], 50, order=8))    #print("90% area distance from coordinate zero: ", get_area_distance(file_location, current_superevent_array, [0,0], 90, order=8))
    cf_list=get_all_null_coords(current_superevent_array[2], which_detectors(current_superevent_array[3]))
    #print("cf_list:", cf_list)
    res_sep_min_50=np.zeros((3, 4))
    res_sep_min_90=np.zeros((3, 4))
    for i in range(3):
        for j in range(4):
            if len(cf_list[i]) !=0:
                if len(cf_list[i][j]) != 0:
                    current_point=cf_list[i][j]
                    res_sep_min_50[i][j]=get_area_distance(file_location, current_superevent_array, current_point, 50, order=8)
                    res_sep_min_90[i][j]=get_area_distance(file_location, current_superevent_array, current_point, 90, order=8)
    #print("angular separations from 90% area: ", res_sep_min_90.tolist())
    sep_90_list.append(res_sep_min_90.tolist())
    sep_50_list.append(res_sep_min_50.tolist())
    #print("angular separations from 50% area: ", res_sep_min_50.tolist())
    #print("Event name: ", current_superevent_array[0])
    #print("RES.OUT")
    iterator0 = iterator0+1
results_d = {'Event name':nameslist, '90% area':area_90_list, '50% area':area_50_list, 'distances from 90% area':sep_90_list, 'distances from 50% area':sep_50_list}
results_df = pd.DataFrame(data=results_d)
print(results_df)
results_df.to_excel("May_events_1_AD_test.xlsx") 
        
