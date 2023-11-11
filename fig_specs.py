import numpy as np
import matplotlib.pyplot as plt
from Calc_Q_fun import calculate_source_spectrum


def fig_specs(freq_vec, channels, station_name, magnitude_p, 
            magnitude_s, mean_mag, event_loc, event_time, 
            seismic_moments_p, seismic_moments_s, 
            average_p_spec, average_s_spec, ajustes,
            deltas):

    n = len(freq_vec)
    if n == 6: 
        print("Plotting both P and S waves spectras")
    else:
        print("Plotting only P or S waves spectras") 

    fig, axes = plt.subplots(1, 4, figsize=(12, 4))
    colors =['blue', 'red', 'green']

    for i in range(n):
        if i <= 2:
            mag_p = np.around(magnitude_p[i], 1)
            axes[i].loglog(freq_vec[i], seismic_moments_p[i], lw=1.0, zorder=2, color=colors[i], label=f'P, Mw={mag_p}')
            axes[i].legend(title= f"Component: {channels[i]}")
            axes[i].grid(True, which='both', zorder=0, alpha=0.5)
            axes[i].set_xlabel("Frequency [Hz]")
            axes[i].set_xlim(np.min(freq_vec[i]), np.max(freq_vec[i]))
            axes[i].set_title(f"Station id: {station_name}")
            if i == 0:
                axes[i].set_ylabel("Seismic moment [$M_O$]")

        if i > 2:
            pos = i - 3
            mag_s = np.around(magnitude_s[pos], 1)
            axes[pos].loglog(freq_vec[i], seismic_moments_s[pos], lw=2, zorder=2, color=colors[pos], label=f'S, Mw={mag_s}')
            axes[pos].legend(title= f"Component: {channels[pos].split('EP')[1]}")
        

    axes[3].loglog(freq_vec[0], average_p_spec, color='k', lw=1, label='mean P', zorder=1)
    axes[3].loglog(freq_vec[0], average_s_spec, color='k', lw=2, label='mean S', zorder=1)
    axes[3].grid(True, which='both', zorder=0, alpha=0.5)
    axes[3].set_title(f"Mean moment rate spectra")
    axes[3].set_xlim(np.min(freq_vec[0]), np.max(freq_vec[0]))
        
    
    for k, ajuste in enumerate(ajustes):
        labels=["$Q_p$", "$Q_s$"]
        Q = np.around(ajuste[0][1], 2)
        cf = np.around(ajuste[2], 2)
        y_model = ajuste[3]
        freq_model = freq_vec[0][:len(y_model)]

        axes[3].loglog(freq_model, y_model, label=f"{labels[k]} = {Q}", zorder=2)
    axes[3].axvline(x=cf, ymin=0, ymax=1, color='magenta', zorder=0, lw=2, label=f"$C_f$= {cf} Hz")
    axes[3].legend()
    axes[3].set_xlabel("Frequency [Hz]")

    fig.suptitle(f'Event_time:{event_time}\nMean Mw = {mean_mag}\nevent_loc = {event_loc}', fontsize=15)

    plt.tight_layout()
    plt.show()
    fig.savefig(f"Spectrum_{station_name}.png", format='png', dpi=600)