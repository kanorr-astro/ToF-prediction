#Utility module for Project Prediction

import math as m
import numpy as np
import matplotlib.pyplot as plt

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
    
    def dist_to_DU(self, dist) ->float: #km to DU
        return dist / 6378.145
    
    def vel_to_DUTU(self, vel) ->float: #km/s to DU/TU
        return vel / 7.90536828
    
    def degs_to_radTU(self, ang_dot) ->float: #deg/s to rad/TU
        return m.radians(ang_dot) / (1.239446309*(10**(-3)))

    def calc_sposition(self) ->list:
        ae = 1                  #Earth Equatorial Radius (DU)
        e = 0.08182             #Earth Oblate Spheroid Eccentricity
        self.xy_dist = abs((ae/m.sqrt(1 - (e**2)*(m.sin(m.radians(self.lat))**2)) + self.dist_to_DU(self.height))*m.cos(m.radians(self.lat)))
        self.z = abs((ae/m.sqrt(1 - (e**2)*(m.sin(m.radians(self.lat))**2)) + self.dist_to_DU(self.height))*m.sin(m.radians(self.lat)))
        return [self.xy_dist*m.cos(m.radians(self.slong)), self.xy_dist*m.sin(m.radians(self.slong)), self.z]
    
    def calc_svel(self) ->list:
        omega_vector = [0, 0, 0.0588336565] #angular momentum vector, rad/TU
        return np.cross(omega_vector, self.sposition)
    
    def track_to_state(self, slant, slant_dot, elev, elev_dot, azim, azim_dot) ->list:
        #Conversions to cononical units and rads
        lat = m.radians(self.lat)
        long = m.radians(self.slong)
        RS = self.sposition
        VS = self.svelocity
        rho = self.dist_to_DU(slant)
        rho_dot = self.vel_to_DUTU(slant_dot)
        el = m.radians(elev)
        el_dot = self.degs_to_radTU(elev_dot)
        az = m.radians(azim)
        az_dot = self.degs_to_radTU(azim_dot)

        #Generate SEZ track vectors and rotation matrix
        rho_vect = ([-rho*m.cos(el)*m.cos(az)],[rho*m.cos(el)*m.sin(az)],[rho*m.sin(el)])
        rho_dot_vect = ([-rho_dot*m.cos(el)*m.cos(az)+rho*m.sin(el)*el_dot*m.cos(az)+rho*m.cos(el)*m.sin(az)*az_dot],
                    [rho_dot*m.cos(el)*m.sin(az)-rho*m.sin(el)*el_dot*m.sin(az)+rho*m.cos(el)*m.cos(az)*az_dot],
                    [rho_dot*m.sin(el)+rho*m.cos(el)*el_dot])
        rotation_matrix = ([m.sin(lat)*m.cos(long), -m.sin(long), m.cos(lat)*m.cos(long)],
                       [m.sin(lat)*m.sin(long), m.cos(long), m.cos(lat)*m.sin(long)],
                       [-m.cos(lat), 0, m.sin(lat)])
        
        #Generate R,V vectors for tracked satellite
        r_sat_ijk = np.dot(rotation_matrix, rho_vect) + ([RS[0]],[RS[1]],[RS[2]])
        v_sat_ijk = np.dot(rotation_matrix, rho_dot_vect) + ([VS[0]],[VS[1]],[VS[2]])
        return [r_sat_ijk, v_sat_ijk]

class Track():
    def __init__(self, state_vector, mu, RE, name:str):
        #Inputs and direct derivatives:
        self.name = name                    #Object name
        self.mu = mu                        #Gravitational Parameter (DU^3 / TU^2)
        self.RE = RE                        #Radius of orbitted body (DU)
        self.R0 = state_vector[0]           #Initial position vector (DU)
        self.V0 = state_vector[1]           #Initial velocity vector (DU/TU)
        self.r0 = self.magnitude(self.R0)   #Initial Radius (DU)
        self.v0 = self.magnitude(self.V0)   #Initial Speed (DU/TU)

        #Orbital Parameters:
        #self.H = np.cross(self.R0, self.V0)
        #self.h = self.magnitude(self.H)
        self.SME = self.spec_mech_energy()

        #Misc:
        self.wireframe = self.wireframe_earth() #[px, py, pz]

    def wireframe_earth(self):
        u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        ae = 1
        px = ae*np.cos(u)*np.sin(v)
        py = ae*np.sin(u)*np.sin(v)
        pz = ae*np.cos(v)
        return [px, py, pz]

    def plot_r0v0(self):
        fig = plt.figure(figsize=(15,15))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim(-1.5,1.5)
        ax.set_ylim(-1.5,1.5)
        ax.set_zlim(-1.5,1.5)
        ax.plot_wireframe(self.wireframe[0], self.wireframe[1], self.wireframe[2], rstride=1, cstride=1, color='green', linewidth=0.5)
        ax.quiver(0,0,0,self.R0[0],self.R0[1],self.R0[2])
        ax.quiver(self.R0[0],self.R0[1],self.R0[2],self.V0[0],self.V0[1],self.V0[2], color='red')
        plt.title('Position and Velocity Vector of Object %s.' %(self.name))
        ax.view_init(azim=0, elev=30)
        plt.show()

    def magnitude(self, vector) ->float:
        return m.sqrt(sum(pow(element, 2) for element in vector))
    
    def spec_mech_energy(self) ->float:
        return ((self.v0**2 / 2) - (self.mu / self.r0))

#    def orb_parameters(self) ->list:
        