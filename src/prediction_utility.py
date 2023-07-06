#Utility module for Project Prediction

import math as m
import numpy as np

class Radar():
    def __init__(self, latitude:float, longitude:float, height:float, date_list:list, time_list:list):
        self.lat = latitude     #Geodetic latitude (degrees)
        self.long = longitude   #Longitude(East +, West -) (degrees)
        self.height = height    #Site altitude above mean sea level (km)
        self.date = date_list   #Date [day, month, year]
        self.time = time_list
        self.UT = time_list[0] + (time_list[1]/60) + (time_list[2]/3600)
        self.JD = (367*self.date[2]) - int((7*(self.date[2] + (int(self.date[1] + 9)/12)))/4) + int((275*self.date[1])/9) + self.date[0] + 1721013.5 + (self.UT/24)
        self.GMST = self.calc_GMST()
        self.slong = self.calc_slong()
        self.sposition = self.calc_sposition()
        self.svelocity = self.calc_svel()
        
    def calc_GMST(self) ->float:
        return ((18.697374558 + 24.06570982441908*(self.JD - 2451545)) % 24)
    
    def calc_slong(self) ->float:
        slong = (((((self.GMST / 24)*360) + self.long) % 360) + 360) % 360
        return slong
    
    def dist_to_DU(self, dist) ->float:
        return dist / 6378.145

    def calc_sposition(self) ->list:
        ae = 1                  #Earth Equatorial Radius (DU)
        e = 0.08182             #Earth Oblate Spheroid Eccentricity
        self.xy_dist = abs((ae/m.sqrt(1 - (e**2)*(m.sin(m.radians(self.lat))**2)) + self.dist_to_DU(self.height))*m.cos(m.radians(self.lat)))
        self.z = abs((ae/m.sqrt(1 - (e**2)*(m.sin(m.radians(self.lat))**2)) + self.dist_to_DU(self.height))*m.sin(m.radians(self.lat)))
        return [self.xy_dist*m.cos(m.radians(self.slong)), self.xy_dist*m.sin(m.radians(self.slong)), self.z]
    
    def calc_svel(self) ->list:
        omega_vector = [0, 0, 0.0588336565] #angular momentum vector, rad/TU
        return np.cross(omega_vector, self.sposition)