import numpy as np
import scipy


def calculate_source_spectrum(frequencies, omega_0, corner_frequency, Q, traveltime, type='RA', gama = 2):
    """
    After Abercrombie (1995) and Boatwright (1980).
    Abercrombie, R. E. (1995). Earthquake locations using single-station deep
    borehole recordings: Implications for microseismicity on the San Andreas
    fault in southern California. Journal of Geophysical Research, 100,
    24003–24013.
    Boatwright, J. (1980). A spectral theory for circular seismic sources,
    simple estimates of source dimension, dynamic stress drop, and radiated
    energy. Bulletin of the Seismological Society of America, 70(1).
    The used formula is:
        Omega(f) = (Omege(0) * e^(-pi * f * T / Q)) / (1 + (f/f_c)^4) ^ 0.5
    :param frequencies: Input array to perform the calculation on.
    :param omega_0: Low frequency amplitude in [meter x second].
    :param corner_frequency: Corner frequency in [Hz].
    :param Q: Quality factor.
    :param traveltime: Traveltime in [s].
    """
    if type == 'RA':
        num = omega_0 * np.exp(-np.pi * frequencies * traveltime / Q)
        denom = (1 + (frequencies / corner_frequency) ** (2*gama)) ** (1/gama)
    elif type == 'Brune':
        num = omega_0
        denom = (1 + (frequencies / corner_frequency) ** (2*gama)) ** (1/gama)
    elif type == 'BW':
        num = omega_0 * np.exp(-np.pi * frequencies * traveltime / Q)
        denom = (1 + (frequencies / corner_frequency) ** 4) ** 0.5

    return num / denom


def get_phase_velocity(pick, V_P,V_S):

    # Only p phase picks.
    if pick.lower() == "p":
        radiation_pattern = 0.52
        velocity = V_P
        k = 0.32
    elif pick.lower() == "s":
        radiation_pattern = 0.63
        velocity = V_S
        k = 0.21
    return radiation_pattern, velocity, k



def fit_Q(spectrum0, frequencies0, traveltime, f_c0, initial_Q, minf, maxf, gama, typeM):
    """
    Fit a theoretical source spectrum to a measured source spectrum.
    Uses a Levenburg-Marquardt algorithm.
    :param spectrum0: The measured source spectrum.
    :param frequencies0: The corresponding frequencies.
    :param traveltime: Event traveltime in [s].
    :param f_c0: the corner frequency
    :param initial_Q: Initial guess for Q.
    :param minf: min frequency for search of Q.
    :param maxf: max frequency for search of Q.
    :param gama: gama valus (use 2)
    :param typeM: is the model type (use 'Brune', 'RA', 'BW' - 'RA' is prefered)
    :returns: Best fits and standard deviations.
    
        popt : [Omega_0, Q]

    """
    omega_1 = np.max(spectrum0[frequencies0 < 20])

    ok = 1
    posO = np.argmin(np.abs(spectrum0 - omega_1))
    f_Omega1 = frequencies0[posO]
    I = np.logical_and(frequencies0 >= minf, frequencies0 <= maxf)
    frequencies = frequencies0[I]
    spectrum = spectrum0[I]

    def f1(frequencies, Q):
            return calculate_source_spectrum(frequencies, omega_1, f_c0, Q, traveltime, typeM, gama)
    
    def f(frequencies, omega_0, Q):
            return calculate_source_spectrum(frequencies, omega_0, f_c0, Q, traveltime, typeM, gama)

    W = np.linspace(100, 1, len(frequencies))
    popt, pcov = scipy.optimize.curve_fit(f1, frequencies, spectrum, sigma=W, p0=initial_Q); popt = [omega_1, popt[0]]  #
    

    if popt is None:
        ok = 0
        return None
    if popt[0] < 0:
        ok = 0
        return None

    if ok == 1:
        Y = calculate_source_spectrum(frequencies, popt[0], f_c0, popt[1], traveltime, typeM, gama)
        y = spectrum
        residuals = np.abs(Y - y)
        residual = np.sqrt(sum(residuals ** 2)) / len(residuals)
        R2 = 1 - np.sum((y-Y)**2) / np.sum((y-np.mean(y))**2)

    return popt, pcov, R2, residuals, Y