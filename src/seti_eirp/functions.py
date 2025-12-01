import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import yaml

# ===============================
# --- Basic Telescope Physics ---
# ===============================

def calc_DishArea(d):
    """Compute dish area in m^2 from diameter in m."""
    return np.pi * (d / 2.0) ** 2


def calc_BeamSize(d, v):
    """
    Compute beam solid angle (approx) in deg^2.

    d : dish diameter [m]
    v : observing frequency [Hz]
    """
    c = 2.998e8  # m/s
    fwhm_deg = 1.22 * (c / (d * v)) * 57.2958
    return (fwhm_deg / 2.0) ** 2 * np.pi


def calc_SEFD(A, Tsys, eff=1.0):
    """
    Calculate SEFD in Jy.

    A    : collecting area [m^2]
    Tsys : system temperature [K]
    eff  : aperture efficiency
    """
    kb = 1.3806488e3  # Boltzmann in Jy·m^2/K
    Ae = A * eff
    return 2.0 * Tsys * kb / Ae


def calc_Sensitivity(m, nu, t, SEFD=0.0, Tsys=10.0, eff=1.0,
                     A=100.0, npol=2.0, narrow=True):
    """
    Radiometer equation (gives sensitivity in Jy).

    m   : S/N threshold
    nu  : channel bandwidth [Hz]
    t   : observing time [s]
    SEFD: system equivalent flux density [Jy] (if 0, computed from A/Tsys/eff)
    """
    if SEFD:
        sefd = SEFD
    else:
        sefd = calc_SEFD(A, Tsys, eff=eff)

    if narrow:
        sens = m * sefd * np.sqrt(nu / (npol * t))
    else:
        sens = m * sefd / np.sqrt(npol * nu * t)

    return sens


def calc_EIRP_min(d_m, Sens):
    """
    Minimum detectable EIRP [W].

    d_m : distance [m]
    Sens: flux density sensitivity [Jy]
    1 Jy = 1e-26 W/m^2/Hz
    """
    return 4.0 * np.pi * d_m**2 * Sens * 1e-26


# ===============================
# --- YAML I/O ---
# ===============================

def load_survey_yaml(path):
    """Load surveys from YAML file."""
    with open(path, "r") as f:
        data = yaml.safe_load(f)
    # Expect top-level key "surveys"
    return data.get("surveys", [])


# ===============================
# --- Survey Calculations ---
# ===============================

def _distance_to_meters(survey):
    """Convert distance fields in YAML to meters."""
    if "max_distance_ly" in survey:
        d_ly = float(survey["max_distance_ly"])
        return d_ly * u.lyr.to("m")
    elif "max_distance_pc" in survey:
        d_pc = float(survey["max_distance_pc"])
        d_ly = d_pc * 3.26156  # pc -> ly
        return d_ly * u.lyr.to("m")
    else:
        raise ValueError(f"No distance (ly or pc) given for survey: {survey.get('name', 'UNKNOWN')}")


def compute_survey_limits(survey):
    """
    Compute EIRP, rarity etc. from one survey dict.

    Expected keys (YAML):
      - name
      - N_stars
      - band_Hz
      - central_freq_Hz
      - dish_diam_m
      - dish_Tsys_K
      - aperture_efficiency
      - SNR_threshold
      - spectral_resolution_Hz
      - obs_time_s
      - max_distance_ly or max_distance_pc
      - instantaneous_bandwidth_Hz
      - SEFD_Jy (optional, if missing computed from dish)
      - npol (optional, default 2)
      - narrow (optional, default True)
    """
    N_stars = float(survey["N_stars"])
    band = float(survey["band_Hz"])
    central_freq = float(survey["central_freq_Hz"])
    dish_diam = float(survey["dish_diam_m"])
    Tsys = float(survey["dish_Tsys_K"])
    eff = float(survey.get("aperture_efficiency", 1.0))
    m = float(survey["SNR_threshold"])
    nu = float(survey["spectral_resolution_Hz"])
    t = float(survey["obs_time_s"])
    iband = float(survey["instantaneous_bandwidth_Hz"])
    npol = float(survey.get("npol", 2.0))
    narrow = bool(survey.get("narrow", True))

    # frequency-normalized rarity term
    freq_range_norm = band / central_freq
    rarity = N_stars * freq_range_norm

    # SEFD
    sefd = survey.get("SEFD_Jy", None)
    if sefd is None:
        A = calc_DishArea(dish_diam)
        sefd = calc_SEFD(A, Tsys, eff=eff)
    else:
        sefd = float(sefd)

    # Sensitivity
    A_for_sens = calc_DishArea(dish_diam)
    Sens = calc_Sensitivity(
        m=m,
        nu=nu,
        t=t,
        SEFD=sefd,
        Tsys=Tsys,
        eff=eff,
        A=A_for_sens,
        npol=npol,
        narrow=narrow,
    )

    # Distance
    d_m = _distance_to_meters(survey)
    EIRP = calc_EIRP_min(d_m, Sens)

    # Sky coverage proxy
    sky = N_stars * calc_BeamSize(dish_diam, central_freq)

    speed = sefd**2 * nu / iband

    return {
        "EIRP": EIRP,
        "rarity": rarity,
        "speed": speed,
        "sky": sky,
        "sefd": sefd,
        "sensitivity": Sens,
    }


# ===============================
# --- Plotting ---
# ===============================
def check_latex():
    """
    Returns True if LaTeX is installed and functional.
    Works for Matplotlib 3.4+ and modern versions without checkdep_usetex.
    """
    import io

    try:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, r"$\alpha+\beta=1$", fontsize=8)
        buf = io.BytesIO()
        fig.savefig(buf, format="pdf")
        plt.close(fig)
        print("LaTeX is available for plotting.... using LaTeX!")
        return True
    except Exception:
        print("LaTeX is not available for plotting.")
        return False

def plot_surveys(yaml_path, output="SETI_limits_comparison.pdf", publish=False, twocolumn=False, noshow=False):
    """
    Load all surveys from YAML, compute EIRP/rarity, and plot them together.
    """
    surveys = load_survey_yaml(yaml_path)
    
   
    if publish:
        import scienceplots
        
        latex = check_latex()
        if latex:
            plt.style.use(['science', 'ieee'])
        else:
            latex = False
            plt.style.use(['science', 'ieee', 'no-latex'])

    
    plt.figure(figsize=(16 if not twocolumn else 9, 9))

    used_names = set()

    for survey in surveys:

        # Handle fixed points
        if survey.get("point_type") == "fixed":
            label = survey["name"] if survey["name"] not in used_names else None
            used_names.add(survey["name"])
            
            if isinstance(label, str) and latex:
                label = (
                    label.replace("\\&", "&")   
                        .replace("&", "\\&")  
                )
            
            edgecolor = survey.get("edgecolor", None)

            plt.plot(
                [survey["logEIRP"]],
                [survey["logRarity"]],
                survey.get("marker", "o"),
                color=survey.get("color", "k"),
                markersize=survey.get("markersize", 16),
                markeredgecolor=edgecolor if edgecolor is not None else None,
                label = label
            )

            continue

        # Computed surveys
        res = compute_survey_limits(survey)
        logE = np.log10(res["EIRP"])
        logR = np.log10(1.0 / res["rarity"])

        # Only label first time name appears
        label = survey["name"] if survey["name"] not in used_names else None
        used_names.add(survey["name"])
        
        if isinstance(label, str) and latex:
            label = (
                label.replace("\\&", "&")   
                    .replace("&", "\\&")  
                )

        plt.plot(
            [logE],
            [logR],
            survey.get("marker", "o"),
            color=survey.get("color", "k"),
            markersize=survey.get("markersize", 16),
            markeredgecolor="w",
            label=label
        )
        
    # -------------------------------------------------------
    # ET Power-law (grey slope)
    # -------------------------------------------------------
    P = np.logspace(10, 26, 100)
    alpha = 0.7 
    
    E1 = 5e11
    S1 = 350
    E2 = 1.98792219e21
    S2 = 7.14285714e8
    alpha = np.log10(S2/S1) / np.log10(E2/E1)
    No = E1**alpha / S1
    NP = No * (1.0 / P**alpha)

    plt.plot(
        np.log10(P),
        np.log10(NP),
        color='gray',
        linewidth=8,
        alpha=0.3
    )
    
    # -------------------------------------------------------
    # Arecibo planetary radar + K1 lines
    # -------------------------------------------------------
    plt.plot([17, 17], [-13, 6], '--', lw=4, color='black', alpha=0.6)
    plt.plot([13, 13], [-13, 6], '-.', lw=4, color='black', alpha=0.6)

    plt.xlabel(r"EIRP $[\log_{10}(\mathrm{Watts})]$", fontsize=20)
    plt.ylabel(
        r"Transmitter Rate "
        r"$\left[ \log_{10} \left( \frac{1}{N_{\mathrm{stars}} \cdot \nu_{\mathrm{rel}}} \right) \right]$"
    , fontsize=20)
    plt.xlim(10, 25)
    plt.ylim(-13, 5)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    # plt.legend(loc=1, prop={"size": 14 if not twocolumn else 6}, labelspacing=1, frameon=True)
    
    # Expand legend across full width and increase spacing
    legend = plt.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.10),
        ncol=6,
        frameon=False,
        labelspacing=1,
        fontsize=12,
        columnspacing=2.0,
        handlelength=2.0,
    )

    legend._legend_box.align = "center"


    plt.subplots_adjust(bottom=0.28) 

    plt.tight_layout()
    plt.savefig(output, format="pdf", dpi=300)
    
    if not noshow:
        plt.show()
        
    print(f"Saved plot to {output}")
