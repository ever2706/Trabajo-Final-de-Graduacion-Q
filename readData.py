import os
from obspy.core import read, UTCDateTime
from processing_st import processing_st
import numpy as np 
from fig_specs import fig_specs
from Calc_Q_fun import fit_Q



def read_nodal_data(data_path, station_id, station_file, event_loc, phase_type, p_arrival, s_arrival, event_time):

    data_full_path = os.path.join(data_path, f"{station_id}..0.*miniseed")

    if phase_type == "P":

        ini_time = UTCDateTime(p_arrival)
        t_before = 0.5
        t_after = 3.0
        st = read(data_full_path, starttime=ini_time-t_before, endtime=ini_time+20)
        st.plot(method='full')
        #Freqs, Specs = processing_st(st, station_file, station_id)

    elif phase_type == "S":

        t_before = 0.5
        t_after = 3.0
        ini_time = UTCDateTime(s_arrival)
        st = read(data_full_path, starttime=ini_time-t_before, endtime=ini_time+t_after)
        #Freqs, Specs = processing_st(st, station_file, station_id)
    
    elif phase_type == "P&S":

        t_before = 0.5
        t_after = 3.0
        ini_time1 = UTCDateTime(p_arrival)
        ini_time2 = UTCDateTime(s_arrival)
        st1 = read(data_full_path, starttime=ini_time1-t_before, endtime=ini_time1+t_after)
        st2 = read(data_full_path, starttime=ini_time2-t_before, endtime=ini_time2+t_after)
        

        # Computing the Spectrums and frequency vectors
        print(f"Processing data for the P-waves now at station {station_id}")
        Freqs1, Channels1, Disp_Specs1, Magnitudes_p, seismic_moments_spec1 = processing_st(
                                                                                        st1, 
                                                                                        station_file=station_file, 
                                                                                        station_id=station_id, 
                                                                                        event_loc=event_loc
                                                                                        )
        print(f"Processing data for the s-waves now at station {station_id}")
        Freqs2, Channels2, Disp_Specs2, Magnitudes_s, seismic_moments_spec2 = processing_st(
                                                                                        st2, 
                                                                                        station_file=station_file, 
                                                                                        station_id=station_id, 
                                                                                        event_loc=event_loc
                                                                                        )

        # Unifying the spectrums and frequency vectors
        total_freqs = Freqs1 + Freqs2
        total_disp_specs = Disp_Specs1 + Disp_Specs2 
        
        # Averaging magnitudes for P and S
        Total_mags = Magnitudes_p + Magnitudes_s
        Total_mags = np.array(Total_mags)
        mean_mag = np.around(np.mean(Total_mags), 1)
        std_mag = np.around(np.std(Total_mags), 1)
        print(f"Mean magnitude computed for this station: {mean_mag} +/- {std_mag}")

        # Averaging the spectrums and seismic moment spectrums
        Mspec_p = np.mean(Disp_Specs1, axis=0)
        Mmoment_p = np.mean(seismic_moments_spec1, axis=0)
        Mspec_p = np.mean(Disp_Specs1, axis=0)

        Mspec_s = np.mean(Disp_Specs2, axis=0)
        Mmoment_s = np.mean(seismic_moments_spec2, axis=0)
        Mspec_s = np.mean(Disp_Specs2, axis=0)
        
        delta_p = ini_time1 - UTCDateTime(event_time)
        delta_s = ini_time2 - UTCDateTime(event_time)

        deltas = [delta_p, delta_s]
        #spectrums = [Mspec_p, Mspec_s]
        spectrums = [Mmoment_p, Mmoment_s]

        #Fitting Q for P wave spectras first
        Gfit = []
        for j, d in enumerate(deltas):
            ajustes =[]
            rs = []
            corners =[]
            synthethics =[]
            for c in np.arange(5,7,0.05):
                popt, pcov, R2, residuals, Y = fit_Q(
                    spectrum0=spectrums[j], 
                    frequencies0=Freqs1[0], 
                    traveltime= deltas[j], 
                    f_c0=c, 
                    initial_Q=300, 
                    minf=0.9, 
                    maxf=25, 
                    gama=2, 
                    typeM='RA')
                
                ajustes.append(popt)
                rs.append(R2)
                corners.append(c)
                synthethics.append(Y)
            
            arg = np.argmax(rs)
            mejor_ajuste = ajustes[arg]
            mejor_r = rs[arg]
            fc = corners[arg]
            model = synthethics[arg]
            

            fitting_vals = [mejor_ajuste, mejor_r, fc, model]
            Gfit.append(fitting_vals)
            
            
        fig_specs(
            freq_vec=total_freqs, 
            channels=Channels1, 
            station_name=station_id, 
            magnitude_p=Magnitudes_p, 
            magnitude_s=Magnitudes_s,
            mean_mag=mean_mag, 
            event_loc=event_loc, 
            event_time=event_time, 
            seismic_moments_p=seismic_moments_spec1,
            seismic_moments_s=seismic_moments_spec2,
            average_p_spec=Mmoment_p, 
            average_s_spec=Mmoment_s,
            ajustes=Gfit, 
            deltas=deltas)

        

