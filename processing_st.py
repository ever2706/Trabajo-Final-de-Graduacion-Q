import numpy as np 
from obspy.core import read, UTCDateTime
from multitaper import MTSpec
from compute_seismic_moment import compute_seismic_moment
from utils import create_response, event_station_distance


def processing_st(st, station_file, station_id, event_loc):

    response = create_response(sensor_model='IGU-16HR-3C')

    Freqs = []
    Specs = []
    Disp_Specs = []
    Channels = []
    Magnitudes = []
    seismic_moments_spec = []
    for i, tr in enumerate(st):
        dt = tr.stats.delta
        channel = tr.stats.channel
        nw = 4.0
        kspec = 7
        tr.detrend('demean')
        tr.detrend('linear')
        tr.simulate(paz_remove=response)
        tr.data *= 100
        tr.filter('bandpass', freqmin=0.5, freqmax=45)
        amp = tr.data
        P = MTSpec(amp, nw, kspec, dt)
        freq, spec = P.rspec()

        sel_freq_range = np.where(np.logical_and(freq >=0.6, freq <= 25))
        freq = freq[sel_freq_range]
        spec = spec[sel_freq_range]

        # Computing the displacement spectrum (dividing by 2*pi*f)
        omega = 2 * np.pi * freq
        disp_spec = spec / omega

        dist = event_station_distance(event_loc, station_file, station_id)

        shear_vel = 3000
        density=2800

        seismic_moment, moment_magnitude = compute_seismic_moment(
                                displacement_spectrum=disp_spec, 
                                frequency=freq, 
                                fmin=1, 
                                fmax=5, 
                                density=density, 
                                shear_velocity=shear_vel, 
                                distance=dist)
        
        seismic_moments = 4 * np.pi * density * shear_vel**3 * dist * np.abs(disp_spec)

        Freqs.append(freq)
        Disp_Specs.append(disp_spec)
        Channels.append(channel)
        Magnitudes.append(moment_magnitude)
        seismic_moments_spec.append(seismic_moments)

    return Freqs, Channels, Disp_Specs, Magnitudes, seismic_moments_spec