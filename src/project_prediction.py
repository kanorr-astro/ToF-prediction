#Project Prediction

#Problem Statement:
#Given: Initial position (r0) and initial velocity (v0) (may be provided from radar observations)
#Find:
    #a. Trajectory type/Conic type
    #b. Position (r) and velocity (v) at impact or closest approach (to Earth)
    #c. Time to impact/closest approach
    #d. Change in true anomaly from t0 to impact/close approach (long way or short way)



#Psuedo
#1. Input type (radar obs or r0/v0)
    #1a. Need module to convert radar obs to r0/v0 (simplify module from site-track)
#2. Process r0/v0 to find orbital parameters and initial true anomaly (nu0)
#3. Determine true anomaly for impact or closest approach
#4. Calculate ToF (delT) and (delnu) to impact/closest approach and position vector at that point.

#Additional Task: Create a class to contain functions related to an object being tracked

#Module imports
import prediction_utility as util
import math as m
import numpy as np

#Inputs
lat = 39.007                        #geodetic latitude (degrees)
long = -104.883                  #longitude in east direction(degrees) (Long w is negative)
height = 2.188464                        #altitude above mean sea level (km)
date_list = [1, 1, 1970]            #Date [day, month, year]
time_list = [0, 0, 0]                     #Time [hour, minutes, seconds]
           
#Input - Radar Contact - Skip if inputting r0/v0 directly
radar = util.Radar(lat, long, height, date_list, time_list)

#Input - r0/v0 directly
print("Date/Time: %s / %s / %s ; %s hours" %(date_list[1], date_list[0], date_list[2], radar.UT))
print("Calculated GMST (hours): ", radar.GMST)
print("Calculated GMST (degrees): ", ((radar.GMST / 24)*360))
print("Radar Site Longitude: ", radar.slong)
print("Radar Site Latitude: ", radar.lat)
print("Radar Site Position Vector (DU): ", radar.sposition)
print("Radar Site Velocity Vector (DU/TU): ", radar.svelocity)



