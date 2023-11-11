import numpy as np 

def compute_seismic_moment(displacement_spectrum, frequency, fmin, fmax, density=2800, shear_velocity=3500, distance=1):
    """
    Compute the seismic moment from the displacement spectrum.

    :param displacement_spectrum: The displacement spectrum (numpy array).
    :param frequency: The frequency values corresponding to the displacement spectrum (numpy array).
    :param density: The density of the medium in kg/m³ (default is 2800 kg/m³).
    :param shear_velocity: The shear wave velocity in m/s (default is 3500 m/s).
    :param distance: The distance from the source to the station in meters (default is 1 m).
    :return: The seismic moment in Newton meters (Nm).
    """
    # Convert shear velocity from km/s to m/s
    shear_velocity *= 1000  # now in m/s

    # Find the low-frequency level of the spectrum (Omega_0)
    # Assuming Omega_0 is the mean of the first few low-frequency spectral values
    # We avoid the very first (DC component) by starting from index 1
    sel_flat_level = np.where(np.logical_and(frequency >= fmin, frequency <=fmax))
    frequency = frequency[sel_flat_level]
    displacement_spectrum = displacement_spectrum[sel_flat_level]
    n = len(displacement_spectrum)
    omega_0 = np.mean(np.abs(displacement_spectrum[1:n]))
    
    # Calculate the seismic moment using the formula
    seismic_moment = 4 * np.pi * density * shear_velocity**3 * distance * omega_0
    moment_magnitude = (2 / 3) * (np.log10(seismic_moment) - 9.1)

    return seismic_moment, moment_magnitude