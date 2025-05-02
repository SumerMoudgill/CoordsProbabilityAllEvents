from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from hpmoc import PartialUniqSkymap
from hpmoc.plot import get_wcs, plot, gridplot
import sys
from CoordsProbabilityCFolder.DetectorMinima import *
from CoordsProbabilityCFolder.RankPoints import *
from CoordsProbabilityCFolder.O4aSuperevents import *
def match_coords(ra, dec, z_ra, z_dec, sigfigs, sigfig_final):
    #print(str(sigfigs)+" significant figures")
    n_coords=len(ra)
    current_i=0
    ra_rounded=np.zeros(n_coords)
    dec_rounded=np.zeros(n_coords)
    z_ra_rounded=round(z_ra, sigfigs)
    z_dec_rounded=round(z_dec, sigfigs)
    for i in range(n_coords):
        ra_rounded[i]=round(ra[i].value, sigfigs)
        dec_rounded[i]=round(dec[i].value, sigfigs)
    target_found=False
    all_coords=False
    final_sigfig=False
    if sigfigs==sigfig_final:
        final_sigfig=True
    #print("finding targets")
    while target_found==False and all_coords==False:
        #print(current_i)
        if current_i>n_coords-1:
            all_coords=True
        elif dec_rounded[current_i]==z_dec_rounded and ra_rounded[current_i]==z_ra_rounded:
            #print("target found")
            #print(ra_rounded[current_i]==z_ra_rounded)
            #print("z_ra_rounded: "+str(z_ra_rounded)+" ; z_dec_rounded: "+str(z_dec_rounded)+" ; ra_rounded[i]: "+str(ra_rounded[current_i])+" ; dec_rounded[i]: "+str(dec_rounded[current_i]))
            target_found=True
        else:
            current_i=current_i+1
    #print("ending target finding loop")
    #print("current_i: "+str(current_i))    
    if target_found==True:
        #print("success")
        #print("z_ra_rounded: "+str(z_ra_rounded)+" ; z_dec_rounded: "+str(z_dec_rounded)+" ; ra_rounded[i]: "+str(ra_rounded[current_i])+" ; dec_rounded[i]: "+str(dec_rounded[current_i]))
        return current_i, sigfigs
    else:
        if final_sigfig==True:
            return None, None
            #print("failure")
        else:
            return match_coords(ra, dec, z_ra, z_dec, sigfigs=sigfigs-1, sigfig_final=sigfig_final)


def get_probability_at_coords(skymap, z_coords, coords_sig_figs=3, final_sig_fig=0):
    ra, dec = skymap.coords()
    #print(z_coords)
    target_dec=z_coords[0]
    target_ra=z_coords[1]
    target_i, C_sigfigs=match_coords(ra, dec, target_ra, target_dec, coords_sig_figs, final_sig_fig)
    if target_i==None:
        return 0, 0, 0
    else:
        coord_string=skymap[target_i]
        probability = single_pixel_to_probability(skymap[target_i])
        return probability, C_sigfigs, coord_string

def get_probabilities_at_minima(event_skymap, timeJulian, detectors):
    #print("b")
    probabilities_list=np.zeros((len(detectors), 4)).tolist()
    C_sigfigs_list = np.zeros((len(detectors), 4)).tolist()
    coord_strings_list=np.zeros((len(detectors), 4)).tolist()
    #timeJulian=float(timeJulian)
    current_detector = 0
    flags_list_H1=[]
    flags_list_L1=[]
    flags_list_V1=[]
    if "H1" in detectors:
        flags_list_H1=[[],[],[],[]]
        detector_stats=H1_detector_stats
        c_array=getMinimaAtJDTArray(detector_stats, timeJulian)
        for i in range(4):
            current_coords=c_array[i]
            probabilities_list[current_detector][i], C_sigfigs_list[current_detector][i], coord_strings_list[current_detector][i]=get_probability_at_coords(event_skymap, current_coords)
            if probabilities_list[current_detector][i] != 0:
                flags_list_H1[i]=current_coords.tolist()
        current_detector=current_detector+1
        #print("H1 done")
    if "L1" in detectors:
        flags_list_L1=[[],[],[],[]]
        detector_stats=L1_detector_stats
        c_array=getMinimaAtJDTArray(detector_stats, timeJulian)
        for i in range(4):
            current_coords=c_array[i]
            probabilities_list[current_detector][i], C_sigfigs_list[current_detector][i], coord_strings_list[current_detector][i]=get_probability_at_coords(event_skymap, current_coords)
            if probabilities_list[current_detector][i] != 0:
                flags_list_L1[i]=current_coords.tolist()
        current_detector=current_detector+1
        #print("L1 done")
    if "V1" in detectors:
        flags_list_V1=[[],[],[],[]]
        detector_stats=V1_detector_stats
        c_array=getMinimaAtJDTArray(detector_stats, timeJulian)
        for i in range(4):
            current_coords=c_array[i]
            probabilities_list[current_detector][i], C_sigfigs_list[current_detector][i], coord_strings_list[current_detector][i]=get_probability_at_coords(event_skymap, current_coords)
            if probabilities_list[current_detector][i] != 0:
                flags_list_V1[i]=current_coords.tolist()
        current_detector=current_detector+1
    flags_list = [flags_list_H1, flags_list_L1, flags_list_V1]
        #print("V1 done")
    #if "K1" in detectors:
    #    detector_stats=K1_detector_stats
    #    c_array=getMinimaAtJDTArray(detector_stats, timeJulian)
    #    for i in range(4):
    #        current_coords=c_array[i]
    #        probabilities_list[current_detector][i], C_sigfigs_list[current_detector][i]=get_probability_at_coords(event_skymap, current_coords)
    #    current_detector=current_detector+1
    #print("c")
    return probabilities_list, C_sigfigs_list, coord_strings_list, flags_list

def superevent_to_probabilities(event_array, file_location="./O4a_fits/"):
    superevent_fits=file_location+event_array[0]+"_bayestar_fits_gz_"+str(event_array[1])+".fits"
    superevent_skymap=PartialUniqSkymap.read(superevent_fits, strategy="ligo")
    detectors_list=which_detectors(event_array[3])
    event_time=event_array[2]
    event_time=float(event_array[2])
    #print("event time")
    #print(event_time)
    #print(type(event_time))
    probabilities_list, C_sigfigs_list, C_coords_list, C_flags=get_probabilities_at_minima(superevent_skymap, event_time, detectors_list)
    return probabilities_list, C_sigfigs_list, C_coords_list, C_flags

def superevent_array_map(event_array, event_name, file_location="./O4a_fits/"):
    superevent_fits=file_location+event_array[0]+"_bayestar_fits_gz_"+str(event_array[1])+".fits"
    print("superevent_fits:", superevent_fits)
    #superevent_fits=file_location+event_array[0]+"_fits,"+str(event_array[1])
    superevent_skymap=PartialUniqSkymap.read(superevent_fits, strategy="ligo")
    detectors=which_detectors(event_array[3])
    event_time=float(event_array[2])
    array_map_list=[]
    if "H1" in detectors:
        filename_d=event_name+"_H1_minima"
        #print(filename_d)
        detector_stats=H1_detector_stats
        c_array=getMinimaAtJDTArray(detector_stats, event_time)
        #print(c_array)
        current_coords=c_array
        plot_tuple=CoordsScatter(c_array, superevent_skymap, filename_d)
        #current_detector=current_detector+1
    if "L1" in detectors:
        filename_d=event_name+"_L1_minima"
        #print(filename_d)
        detector_stats=L1_detector_stats
        c_array=getMinimaAtJDTArray(detector_stats, event_time)
        current_coords=c_array
        plot_tuple=CoordsScatter(c_array, superevent_skymap, filename_d)
        #current_detector=current_detector+1     
    if "V1" in detectors:
        filename_d=event_name+"_V1_minima"
        #print(filename_d)
        detector_stats=V1_detector_stats
        c_array=getMinimaAtJDTArray(detector_stats, event_time)
        current_coords=c_array
        plot_tuple=CoordsScatter(c_array, superevent_skymap, filename_d)
        #current_detector=current_detector+1

