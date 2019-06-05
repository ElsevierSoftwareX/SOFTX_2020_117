import numpy as np


def filter_cavity_parameters(P, lambda0, m, Titm, Tsrm, L, loss=0):
    """
    Compute the optimal filter cavity detuning phase and input mirror
    transmissivity for frequency-dependent squeezing in a tuned, dual-recycled
    interferometer.

    P - Arm cavity circulating power [W]
    lambda0 - main laser wavelength [m]
    m - mass of test masses (mirrors) [kg]
    Titm - Power transmissivity of input test mass
    Tsrm - Power transmissivity of signal recycling mirror
    L - Filter cavity length [m]
    loss - Filter cavity round-trip power loss

    See Kwee et al. `Decoherence and degradation of squeezed states in
    quantum filter cavities` Physical Review D, 2014, 90, 062006
    """
    c = 299792458.0
    w0 = 2 * np.pi * c / lambda0
    FSR = c / (2.0 * L)

    # Approximate frequency at which IFO reaches SQL
    # (i.e. radiation pressure noise = shot noise)
    w_sql_0 = 8 / c * np.sqrt(P * w0 / (m * Titm))
    w_sql = np.sqrt(Tsrm) / (1 + np.sqrt(1 - Tsrm)) * w_sql_0

    if loss > 0:
        epsilon = 4 / (
            2 + np.sqrt(2 + 2 * np.sqrt(1 + (2 * w_sql / (FSR * loss)) ** 4))
        )
    else:
        epsilon = 0

    # Filter cavity bandwidth
    gamma_f = (
        np.sqrt(2 / ((2 - epsilon) * np.sqrt(1 - epsilon)))
        * w_sql
        / (2 * np.pi * np.sqrt(2))
    )

    # Detuning
    delta = np.sqrt(1 - epsilon) * gamma_f

    phi = -delta / FSR * 180.0
    T = 2.0 * gamma_f / FSR * 2 * np.pi
    return phi, T
