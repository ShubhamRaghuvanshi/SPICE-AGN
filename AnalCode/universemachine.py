import numpy as np

# Behroozi+2019 / UniverseMachine (arXiv:1806.07893), Appendix J / Table J1
# Row: "Obs., All, All, Excl."
UM_DR1_OBS_ALL_ALL_EXCL = dict(
    epsilon0=-1.435, epsilona=-0.075, epsilon_lna=1.831, epsilonz=-2.818,
    M0=12.035, Ma=-0.100, M_lna=4.556, Mz=-2.453,
    alpha0=4.417, alphaa=-2.255, alpha_lna=-0.731, alphaz=-0.064,
    beta0=1.963, betaa=-0.010, betaz=-2.316,
    delta0=0.482,
    gamma0=-0.841, gammaa=-0.471, gammaz=-1.055,
)

def mstar_from_mpeak_universemachine(Mpeak, z):
    """
    Median stellar mass M* from UniverseMachine (Behroozi+2019).

    Parameters
    ----------
    Mpeak : float or array_like
        Peak halo mass in Msun.
    z : float
        Redshift (scalar).

    Returns
    -------
    Mstar : np.ndarray or float
        Median stellar mass in Msun.
    """
    p = UM_DR1_OBS_ALL_ALL_EXCL

    Mpeak = np.asarray(Mpeak, dtype=float)

    # enforce scalar redshift
    z = float(z)

    a = 1.0 / (1.0 + z)
    lna = np.log(a)

    # log10(M1)
    log10_M1 = (
        p["M0"]
        + p["Ma"] * (a - 1.0)
        - p["M_lna"] * lna
        + p["Mz"] * z
    )
    M1 = 10.0 ** log10_M1

    # epsilon
    eps = (
        p["epsilon0"]
        + p["epsilona"] * (a - 1.0)
        - p["epsilon_lna"] * lna
        + p["epsilonz"] * z
    )

    # alpha
    alpha = (
        p["alpha0"]
        + p["alphaa"] * (a - 1.0)
        - p["alpha_lna"] * lna
        + p["alphaz"] * z
    )

    # beta
    beta = p["beta0"] + p["betaa"] * (a - 1.0) + p["betaz"] * z

    # delta
    delta = p["delta0"]

    # gamma
    log10_gamma = p["gamma0"] + p["gammaa"] * (a - 1.0) + p["gammaz"] * z
    gamma = 10.0 ** log10_gamma

    # x = log10(Mpeak/M1)
    x = np.log10(Mpeak / M1)

    # Eq. (J1)
    log10_Mstar_over_M1 = (
        eps
        - np.log10(10.0 ** (-alpha * x) + 10.0 ** (-beta * x))
        + gamma * np.exp(-0.5 * (x / delta) ** 2)
    )

    return M1 * (10.0 ** log10_Mstar_over_M1)
