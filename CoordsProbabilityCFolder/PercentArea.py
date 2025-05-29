import numpy as np
#import healpy as hp #this does skymap plotting (based on matplotlib)
#from healpy import nside2pixarea
import pandas as pd
import healpy as hp
from astropy.io import fits
from astropy.table import Column, Table, join

from ligo.skymap.moc import rasterize #this processes skymaps 
from ligo.skymap.io.fits import read_sky_map #this reads skymaps

#for i in range(len(O4a_superevents)):
#def get_long_lat_deg(healpy_array, order=8):
#    lat_long_rad_array = [hp.healpix_to_lonlat(healpy_array, 2**order)
#    long_lat_array=np.zeros((len(healpy_array), 2))
#    for i in range(len(healpy_array)):
#        long_lat_array[i][0]=np.radtodeg(long_lat_array[i][1])
#        long_lat_array[i][1]=np.radtodeg(long_lat_array[i][0])
#    return long_lat_array
def get_lat_lon_deg(healpy_array, order=8):
    lat_lon_array=np.zeros((len(healpy_array), 2))
    for i in range(len(healpy_array)):
        current_lon, current_lat = hp.pixelfunc.pix2ang(2**order, i, nest=True, lonlat=True)
        lat_lon_array[i][1]=current_lon
        lat_lon_array[i][0]=current_lat
    return lat_lon_array


def get_nested_data(file_location, event_array, order=8):
    superevent_fits=file_location+event_array[0]+"_bayestar_fits_gz_"+str(event_array[1])+".fits"
    M1_data=read_sky_map(superevent_fits, nest=False, distances=True, moc=True)
    M1_nested_data=rasterize(M1_data,order=order)
        #when rasterize: nside = 2^order 
        #number of pixels = 12*(n_side^2) 
        #order shouldn't go more than 9 since we want it fast. 
    #onepixelarea=nside2pixarea(256, degrees = True)
    #print("nested data length: ", len(M1_nested_data))
    #for i in range(len(M1_nested_data)):
        #print("M1 nested data", i, "is", M1_nested_data[i])
    #M1_probdensity = M1_nested_data ['PROBDENSITY']
    #M1_mu=M1_nested_data['DISTMU']
    #M1_sigma=M1_nested_data['DISTSIGMA']
    #M1_sort_prob=np.sort(M1_probdensity/sum(M1_probdensity))[::-1]
    #print(M1_nested_data)
    return M1_nested_data

def rank_nested_data(file_location, event_array, order=8):
    M1_nested_data=get_nested_data(file_location, event_array, order=order)
    #print("M1_nested_data: ", M1_nested_data)
    #print("M1_nested_data.columns: ", M1_nested_data.columns)
    M1_probdensity = M1_nested_data ['PROBDENSITY']
    M1_mu=M1_nested_data['DISTMU']
    M1_sigma=M1_nested_data['DISTSIGMA']
    M1_sort_prob=np.sort(M1_probdensity/sum(M1_probdensity))[::-1]
    return M1_sort_prob

def rank_combined_data(file_location, event_array, order=8):
    M1_nested_data=get_nested_data(file_location, event_array, order=order)
    M1_coordinates=get_lat_lon_deg(M1_nested_data, order=order)
    if len(M1_coordinates)==len(M1_nested_data):
        M1_combined_data=np.zeros((len(M1_nested_data), 3))
        for i in range(len(M1_nested_data)):
            M1_combined_data[i][0]=M1_nested_data[i][0]
            M1_combined_data[i][1]=M1_coordinates[i][0] #declination
            M1_combined_data[i][2]=M1_coordinates[i][1] #right ascension
        #print("len(M1_combined_data):", len(M1_combined_data))
        #print("M1_combined_data: ", M1_combined_data)
        M1_combined_data_ranked=M1_combined_data[M1_combined_data[:,0].argsort()]
        M1_cdr_flip = np.flip(M1_combined_data_ranked, axis=0)
        return M1_cdr_flip
    else:
        print("lengths do not match")
        print(len(M1_nested_data))
        print(len(M1_coordinates))
    

#def angular_separation_deg_hp(lat_lon_o, point_test):
#    dec_rad=np.deg2rad(point_test[0])
#    ra_rad=np.deg2rad(point_test[1])
#    point_dec_rad=np.deg2rad(lat_lon_o[1])
#    point_ra_rad=np.deg2rad(lat_lon_o[2])
#    sep_rad=angular_separation(dec_rad, ra_rad, point_dec_rad, point_ra_rad)
#    return np.rad2deg(sep_rad)

def angular_separation_deg_hp(lat_lon_o, point_test):
    dec_rad=np.deg2rad(point_test[0])
    ra_rad=np.deg2rad(point_test[1])
    point_dec_rad=np.deg2rad(lat_lon_o[1])
    point_ra_rad=np.deg2rad(lat_lon_o[2])
    sep_rad=np.arccos((np.sin(dec_rad)*np.sin(point_dec_rad))+(np.cos(dec_rad)*np.cos(point_dec_rad)*np.cos(point_ra_rad-ra_rad)))
    sep_deg=np.rad2deg(sep_rad)
    return sep_deg

def angular_separation_array(combined_dataset, test_point, order=8, debug=False):
    #print("combinded_dataset:", combined_dataset)
    separations_list=[]
    for i in range(len(combined_dataset)):
        separations_list.append(angular_separation_deg_hp(combined_dataset[i], test_point))
    #print("separations_list: ", separations_list)
    return separations_list

def get_area_sum(cdr):
    area_sum=0
    for i in range(len(cdr)):
        area_sum=area_sum+cdr[i][0]
    return area_sum

def get_area_distance(file_location, event_array, test_point, area_percent, order=8, debug=False):
    M1_cdr=rank_combined_data(file_location, event_array, order=order)
    M1_count=0
    M1_total=0
    area_value=area_percent/100
    area_sum=get_area_sum(M1_cdr)
    #print("M1_cdr:", M1_cdr)
    for i in range(0,len(M1_cdr)):
        currentval = M1_cdr[i][0]/area_sum
        M1_total=M1_total+currentval
        #print("i:", i)
        #print("M1_total:", M1_total)
        #print("M1_count:", M1_count)
        #print("area_value:", area_value)
        #print("M1_total<=area_value", M1_total<=area_value)
        if M1_total<=area_value:
            M1_count+=1
        else:
            #print("breaking loop")
            break
    #print("final M1_count:", M1_count)
    M1_cdr_reduced=M1_cdr[0:M1_count+1]
    #print("len(M1_cdr_reduced):", len(M1_cdr_reduced))
    angular_sep_reduced=angular_separation_array(M1_cdr_reduced, test_point, order=order)
    print("len(angular_sep_reduced):", len(angular_sep_reduced))
    #print("angular_sep_reduced: ", angular_sep_reduced)
    min_distance=min(angular_sep_reduced)
    return min_distance


def get_area(file_location, event_array, area_percent, onepixelarea=0.052455852825697924):
    M1_sort_prob=rank_nested_data(file_location, event_array)
    M1_count=0
    M1_total=0
    area_value=area_percent/100
    for i in range(0,len(M1_sort_prob)): 
        M1_total+=M1_sort_prob[i]
        if M1_total<=area_value:
            M1_count+=1
        else: 
            break
    return M1_count*onepixelarea
