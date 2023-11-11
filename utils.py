import numpy as np 
from obspy.geodetics.base import gps2dist_azimuth 

def create_response(sensor_model='IGU-16HR-3C'):

    zeros = [14164 + 0.0j, -7162 + 0.0j, 0 + 0.0j, 0 + 0.0j]
    poles = [-1720.4 + 0j, -1.2 + 0.9j, -1.2 - 0.9j]
    response={'poles': poles, 'zeros': zeros, 'gain': 3355.4428, 'sensitivity':76.7}

    return response
    

def event_station_distance(event_loc, station_file, station_id):

    staname, stalat, stalon,  = np.loadtxt(station_file, usecols=(0, 2, 3), dtype='str', unpack=True, skiprows=1)
    stalat = stalat.astype(float)
    stalon = stalon.astype(float)

    sel_station_loc = np.nonzero(staname == station_id)
    sla = stalat[sel_station_loc]
    slo = stalon[sel_station_loc]
    staname = staname[sel_station_loc]

    dist, az, baz = gps2dist_azimuth(lat1=event_loc[0], lon1=event_loc[1], lat2=sla[0], lon2=slo[0])
    dist = np.around(dist, 2)

    return dist