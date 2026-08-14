# ============================================================
# Gaia DR3 + galpy Galactic Spiral-Arm Project
# ============================================================
#
# Scientific question:
#
#   How does the stellar kinematic response of a Galactic disc
#   depend on the pattern speed of a rotating spiral perturbation?
#
# The project combines:
#
#   Gaia DR3 observations
#          +
#   fast galpy test-particle simulations
#          +
#   an exploratory pattern-speed response scan
#
# This is an ongoing independent research project. The present version is
# an exploratory test-particle analysis, not a full SBI inference.
# ============================================================

# ============================================================
# 1. IMPORTS
# ============================================================

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

import astropy.units as u
from astropy.coordinates import SkyCoord, Galactocentric

from astroquery.gaia import Gaia

from galpy.orbit import Orbit
from galpy.potential import (
    MWPotential2014,
    SteadyLogSpiralPotential,
    vcirc,
)

def set_publication_style():

    mpl.rcParams.update({

        "font.family": "serif",
        "font.serif": ["STIXGeneral"],

        "font.size": 13,

        "axes.labelsize": 14,
        "axes.titlesize": 15,

        "xtick.labelsize": 12,
        "ytick.labelsize": 12,

        "legend.fontsize": 11,

        # -------------------------
        # Axes
        # -------------------------

        "axes.linewidth": 1.2,

        "xtick.direction": "in",
        "ytick.direction": "in",

        "xtick.top": True,
        "ytick.right": True,

        "xtick.major.size": 6,
        "ytick.major.size": 6,

        "xtick.minor.size": 3,
        "ytick.minor.size": 3,

        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,

        "xtick.minor.width": 0.8,
        "ytick.minor.width": 0.8,

        "xtick.minor.visible": True,
        "ytick.minor.visible": True,

        # -------------------------
        # Figure
        # -------------------------

        "figure.dpi": 150,
        "savefig.dpi": 400,

        "savefig.bbox": "tight",

        # -------------------------
        # General
        # -------------------------

        "axes.grid": False,

        "legend.frameon": False,

        "lines.linewidth": 1.8,

        "lines.markersize": 6,

    })


set_publication_style()

GAIA_FILE = Path("gaia_dr3_local_rvs.csv")

# Number of Gaia stars to request.
N_GAIA = 30000


# ------------------------------------------------------------
# Galactocentric reference frame
# ------------------------------------------------------------

GALCEN_DISTANCE = 8.122       # kpc
Z_SUN = 0.0208                # kpc

# Solar velocity in Galactocentric Cartesian coordinates.
#
# [vx, vy, vz] in km/s
V_SUN = np.array([
    12.9,
    245.6,
    7.78
])


# ------------------------------------------------------------
# galpy model
# ------------------------------------------------------------

# MWPotential2014 uses the conventional galpy scaling:
#
# R0 = 8 kpc
# V0 = 220 km/s

RO_MODEL = GALCEN_DISTANCE               # kpc
VO_MODEL = 220.0              # km/s


# ------------------------------------------------------------
# Spiral-arm model
# ------------------------------------------------------------

# Recent Gaia-based work finds spiral pattern speeds in the
# broad ~10--20 km/s/kpc range, but we deliberately explore
# a somewhat wider range here.

PATTERN_SPEEDS = np.arange(
    8.0,
    28.0,
    2.0
)

# Two-armed logarithmic spiral
N_ARMS = 2

# Pitch angle
PITCH_DEG = 10.0

# Spiral perturbation strength.
#
# This is intentionally modest; this is a toy model.
SPIRAL_A = -0.02

# Spiral phase relative to the Sun-GC line
GAMMA_DEG = 45.0


# ------------------------------------------------------------
# Numerical settings
# ------------------------------------------------------------

# Number of test particles.
N_SIM_STARS = 400

# Total integration time.
T_END = 0.5

# Number of output points.
N_OUTPUT = 51

# Fixed integration step.
#
# This is deliberately coarse for speed, can try with 0.005 and 0.0025.
DT = 0.01            # Gyr


# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

RNG = np.random.default_rng(42)


# ============================================================
# 3. DOWNLOAD GAIA DR3 DATA
# ============================================================

def query_gaia():

    # Deliberately querying a somewhat larger random-index slice
    # and then applying the astrophysical cuts.

    RANDOM_MAX = 20_000_000

    query = f"""
    SELECT TOP {N_GAIA}

        source_id,
        ra,
        dec,

        parallax,
        parallax_error,
        parallax_over_error,

        pmra,
        pmdec,

        radial_velocity,
        radial_velocity_error,

        ruwe,
        visibility_periods_used,
        phot_g_mean_mag,
        random_index

    FROM gaiadr3.gaia_source

    WHERE random_index < {RANDOM_MAX}

      AND parallax > 0.15
      AND parallax_over_error > 5

      AND radial_velocity IS NOT NULL
      AND radial_velocity_error < 5

      AND ruwe < 1.4
      AND visibility_periods_used >= 10
      AND phot_g_mean_mag < 16.5
    """

    print("\nQuerying Gaia DR3...")
    print("Selecting a representative random-index subset...")

    try:

        job = Gaia.launch_job_async(
            query,
            dump_to_file=False
        )

        table = job.get_results()

    except Exception as exc:

        print("\nGaia archive query failed.")
        print(f"Error: {exc}")
        print(
            "\nRetrying with a smaller random-index slice..."
        )

        RANDOM_MAX_SMALL = 5_000_000

        query_small = f"""
        SELECT TOP {N_GAIA}

            source_id,
            ra,
            dec,

            parallax,
            parallax_error,
            parallax_over_error,

            pmra,
            pmdec,

            radial_velocity,
            radial_velocity_error,

            ruwe,
            visibility_periods_used,
            phot_g_mean_mag,
            random_index

        FROM gaiadr3.gaia_source

        WHERE random_index < {RANDOM_MAX_SMALL}

          AND parallax > 0.15
          AND parallax_over_error > 5

          AND radial_velocity IS NOT NULL
          AND radial_velocity_error < 5

          AND ruwe < 1.4
          AND visibility_periods_used >= 10
          AND phot_g_mean_mag < 16.5
        """

        job = Gaia.launch_job_async(
            query_small,
            dump_to_file=False
        )

        table = job.get_results()

    df = table.to_pandas()

    df.to_csv(
        GAIA_FILE,
        index=False
    )

    print(
        f"Downloaded {len(df):,} Gaia DR3 stars."
    )

    print(
        f"Cached catalogue as {GAIA_FILE}"
    )

    return df


# ============================================================
# 4. TRANSFORM GAIA → GALACTOCENTRIC COORDINATES
# ============================================================

def transform_gaia(df):

    """
    Convert Gaia ICRS astrometry into a Galactocentric
    Cartesian frame and then cylindrical coordinates.
    """

    # --------------------------------------------------------
    # Distance
    # --------------------------------------------------------
    #
    # Because the selected stars have parallax S/N > 10,
    # inverse parallax is adequate for this quick project.
    #
    # This is NOT intended as a precision distance estimator.

    distance = (
        1000.0 / df["parallax"].to_numpy()
    ) * u.pc


    # --------------------------------------------------------
    # Gaia ICRS coordinates
    # --------------------------------------------------------

    coords = SkyCoord(

        ra=df["ra"].to_numpy() * u.deg,

        dec=df["dec"].to_numpy() * u.deg,

        distance=distance,

        pm_ra_cosdec=(
            df["pmra"].to_numpy()
            * u.mas
            / u.yr
        ),

        pm_dec=(
            df["pmdec"].to_numpy()
            * u.mas
            / u.yr
        ),

        radial_velocity=(
            df["radial_velocity"].to_numpy()
            * u.km
            / u.s
        ),

        frame="icrs",
    )

    galactic = coords.galactic

    galactic_latitude = (
        galactic.b.to_value(u.deg)
)
    # --------------------------------------------------------
    # Explicit Galactocentric frame
    # --------------------------------------------------------

    gc_frame = Galactocentric(

        galcen_distance=(
            GALCEN_DISTANCE * u.kpc
        ),

        galcen_v_sun=(
            V_SUN * u.km / u.s
        ),

        z_sun=(
            Z_SUN * u.kpc
        ),
    )


    # Transform
    gc = coords.transform_to(
        gc_frame
    )


    # --------------------------------------------------------
    # Cartesian positions
    # --------------------------------------------------------

    x = gc.x.to_value(
        u.kpc
    )

    y = gc.y.to_value(
        u.kpc
    )

    z = gc.z.to_value(
        u.kpc
    )


    # --------------------------------------------------------
    # Cartesian velocities
    # --------------------------------------------------------

    vx = gc.v_x.to_value(
        u.km / u.s
    )

    vy = gc.v_y.to_value(
        u.km / u.s
    )

    vz = gc.v_z.to_value(
        u.km / u.s
    )


    # --------------------------------------------------------
    # Cylindrical coordinates
    # --------------------------------------------------------

    R = np.hypot(
        x,
        y
    )

    phi = np.arctan2(
        y,
        x
    )


    # Radial velocity
    vR = (
        x * vx
        +
        y * vy
    ) / R


    # Astropy uses a right-handed Galactocentric frame in which the Sun
    # lies at negative x.  In that convention the cylindrical v_phi of
    # prograde disk rotation is negative.  For readability below I store
    # a prograde-positive azimuthal velocity.
    vT = -(
        x * vy
        -
        y * vx
    ) / R


    # --------------------------------------------------------
    # Build output dataframe
    # --------------------------------------------------------

    out = df.copy()

    out["distance_kpc"] = (
        distance.to_value(
            u.kpc
        )
    )

    out["x_kpc"] = x
    out["y_kpc"] = y
    out["z_kpc"] = z

    out["vx_kms"] = vx
    out["vy_kms"] = vy
    out["vz_kms"] = vz

    out["R_kpc"] = R
    out["phi_rad"] = phi
    out["galactic_latitude_deg"] = galactic_latitude

    out["vR_kms"] = vR
    out["vT_kms"] = vT


    print("\nGalactocentric extent before cuts:")

    print(
        f"R: {np.nanmin(out['R_kpc']):.2f} "
        f"to {np.nanmax(out['R_kpc']):.2f} kpc"
    )

    print(
        f"x: {np.nanmin(out['x_kpc']):.2f} "
        f"to {np.nanmax(out['x_kpc']):.2f} kpc"
    )

    print(
        f"y: {np.nanmin(out['y_kpc']):.2f} "
        f"to {np.nanmax(out['y_kpc']):.2f} kpc"
    )

    print(
        f"z: {np.nanmin(out['z_kpc']):.2f} "
        f"to {np.nanmax(out['z_kpc']):.2f} kpc"
    )


    return out


# ============================================================
# 5. FIGURE 1 — GAIA RADIAL VELOCITY FIELD
# ============================================================

def plot_gaia_radial_velocity_map(df):

    plot_df = df[np.isfinite(df["x_kpc"]) & np.isfinite(df["y_kpc"]) & np.isfinite(df["vR_kms"])].copy()

    fig, ax = plt.subplots(figsize=(8.2, 6.8))

    hb = ax.hexbin(
        plot_df["x_kpc"],
        plot_df["y_kpc"],
        C=plot_df["vR_kms"],
        reduce_C_function=np.nanmedian,
        gridsize=55,
        mincnt=3,
        cmap="coolwarm",
        vmin=-30,
        vmax=30,
        linewidths=0.0,
    )

    cb = fig.colorbar(hb, ax=ax, pad=0.02)
    cb.set_label(r"Median $v_R$ [km s$^{-1}$]")

    # Astropy's Galactocentric frame places the Sun at negative x.
    ax.scatter(
        -GALCEN_DISTANCE, 0, marker="*", s=130, color="black",
        label="Sun", zorder=10,
    )

    # Robust limits: show the observed local volume clearly.
    xmin, xmax = np.nanpercentile(plot_df["x_kpc"], [0.5, 99.5])
    ymin, ymax = np.nanpercentile(plot_df["y_kpc"], [0.5, 99.5])
    xpad = max(0.35, 0.08 * (xmax - xmin))
    ypad = max(0.35, 0.08 * (ymax - ymin))
    ax.set_xlim(xmin - xpad, xmax + xpad)
    ax.set_ylim(ymin - ypad, ymax + ypad)

    # Direction to the Galactic centre (which lies outside this local view).
    ax.annotate(
        "toward Galactic centre",
        xy=(0.97, 0.08), xytext=(0.72, 0.08),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=dict(arrowstyle="->", lw=1.2),
        ha="center", va="center", fontsize=10,
    )

    ax.set_xlabel("Galactocentric x [kpc]")
    ax.set_ylabel("Galactocentric y [kpc]")
    ax.set_title("Local Gaia DR3 radial-velocity field")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="upper left")

    fig.tight_layout()
    fig.savefig("figure1_gaia_vR_xy.png", dpi=400, bbox_inches="tight")
    plt.show()


# ============================================================
# 6. FIGURE 2 — GAIA VELOCITY DISTRIBUTION
# ============================================================

def plot_gaia_velocity_distribution(df):

    # A thin-disk-like subset makes the local phase-space structure easier to see.
    plot_df = df[
        np.isfinite(df["vR_kms"])
        & np.isfinite(df["vT_kms"])
        & (np.abs(df["z_kpc"]) < 0.8)
    ].copy()

    # Robust plotting limits prevent a small number of halo/outlier velocities
    # from dominating the axes without deleting them from the underlying data.
    vr_lo, vr_hi = np.nanpercentile(plot_df["vR_kms"], [0.5, 99.5])
    vt_lo, vt_hi = np.nanpercentile(plot_df["vT_kms"], [0.5, 99.5])

    fig, ax = plt.subplots(figsize=(8.0, 6.2))
    hb = ax.hexbin(
        plot_df["vR_kms"],
        plot_df["vT_kms"],
        gridsize=60,
        mincnt=2,
        bins="log",
        cmap="viridis",
        extent=(vr_lo, vr_hi, vt_lo, vt_hi),
        linewidths=0.0,
    )

    cb = fig.colorbar(hb, ax=ax, pad=0.02)
    cb.set_label("Number of stars")

    ax.set_xlim(vr_lo, vr_hi)
    ax.set_ylim(vt_lo, vt_hi)
    ax.axvline(0, color="0.5", lw=0.8, ls="--")
    ax.set_xlabel(r"$v_R$ [km s$^{-1}$]")
    ax.set_ylabel(r"Prograde $v_\phi$ [km s$^{-1}$]")
    ax.set_title("Gaia DR3 local Galactocentric velocity distribution")

    fig.tight_layout()
    fig.savefig("figure2_gaia_velocity_distribution.png", dpi=400, bbox_inches="tight")
    plt.show()


# ============================================================
# 7. FIGURE 3 — GAIA RADIAL VELOCITY VS AZIMUTH
# ============================================================

def plot_gaia_azimuthal_profile(df):

    annulus = df[
        (df["R_kpc"] > 7.5)
        & (df["R_kpc"] < 8.5)
        & np.isfinite(df["phi_rad"])
        & np.isfinite(df["vR_kms"])
    ].copy()

    # The local Gaia sample sits around the Sun at phi ~= +/- pi in the
    # Astropy frame.  Re-centre azimuth on the solar direction so the
    # observed wedge is continuous rather than split across the branch cut.
    annulus["dphi"] = (annulus["phi_rad"] - np.pi + np.pi) % (2 * np.pi) - np.pi

    lo, hi = np.nanpercentile(annulus["dphi"], [1, 99])
    bins = np.linspace(lo, hi, 13)
    centres = 0.5 * (bins[:-1] + bins[1:])

    medians, errors, counts = [], [], []
    for a, b in zip(bins[:-1], bins[1:]):
        values = annulus.loc[(annulus["dphi"] >= a) & (annulus["dphi"] < b), "vR_kms"].to_numpy()
        values = values[np.isfinite(values)]
        counts.append(len(values))
        if len(values) >= 15:
            med = np.median(values)
            # Approximate standard error of a median for a well-behaved sample.
            err = 1.253 * np.std(values, ddof=1) / np.sqrt(len(values))
            medians.append(med)
            errors.append(err)
        else:
            medians.append(np.nan)
            errors.append(np.nan)

    medians = np.asarray(medians)
    errors = np.asarray(errors)

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    ax.errorbar(
        centres, medians, yerr=errors, marker="o", capsize=3,
        label="Gaia DR3",
    )
    ax.axhline(0, linestyle="--", linewidth=1, color="0.45")
    ax.set_xlabel(r"Azimuth relative to the solar direction $\Delta\phi$ [rad]")
    ax.set_ylabel(r"Median $v_R$ [km s$^{-1}$]")
    ax.set_title(r"Gaia DR3 radial-velocity trend in the solar annulus ($7.5<R<8.5$ kpc)")
    ax.legend()

    fig.tight_layout()
    fig.savefig("figure3_gaia_azimuthal_vR.png", dpi=400, bbox_inches="tight")
    plt.show()


# ============================================================
# 8. M = 2 FOURIER MEASURE
# ============================================================
PHI_WEDGE_HALFWIDTH = 0.95
def m2_fourier(
    R,
    phi,
    vR,
    rmin=7.5,
    rmax=8.5,
    phi_center=None,
    phi_halfwidth=None,
):
    mask = (

        (R >= rmin)
        &
        (R <= rmax)
        &
        np.isfinite(vR)
        &
        np.isfinite(phi)
    )
    if phi_halfwidth is not None:

        dphi = (
            (phi - phi_center + np.pi)
            % (2 * np.pi)
            - np.pi
        )

        mask &= np.abs(dphi) < phi_halfwidth


    if np.sum(mask) < 12:

        return np.nan, np.nan


    # m=2 Fourier coefficient
    coefficient = np.mean(

        vR[mask]
        *
        np.exp(
            -2j * phi[mask]
        )
    )


    # Real-space sinusoidal amplitude
    amplitude = (
        2.0
        *
        np.abs(
            coefficient
        )
    )


    phase = (
        0.5
        *
        np.angle(
            coefficient
        )
    )


    return amplitude, phase


# ============================================================
# 9. CREATE SPIRAL POTENTIAL
# ============================================================

def make_spiral_potential(
    pattern_speed
):

    # Convert km/s/kpc → galpy frequency units.
    #
    # galpy's natural frequency scale is V0/R0.

    omega_internal = (

        pattern_speed
        /
        (VO_MODEL / RO_MODEL)
    )


    spiral = SteadyLogSpiralPotential(

        amp=1.0,

        omegas=omega_internal,

        A=SPIRAL_A,

        m=N_ARMS,

        p=np.deg2rad(
            PITCH_DEG
        ),

        gamma=np.deg2rad(
            GAMMA_DEG
        ),

        ro=RO_MODEL,

        vo=VO_MODEL,
    )


    return (
        MWPotential2014
        +
        spiral
    )


# ============================================================
# 10. CREATE A COLD MOCK GALACTIC DISC
# ============================================================

def initial_cold_disc(
    n=N_SIM_STARS
):

    # Radial positions
    R_kpc = RNG.uniform(
        7.0,
        9.0,
        n
    )


    # Uniform azimuth
    phi = RNG.uniform(
        -np.pi,
        np.pi,
        n
    )


    # Convert R to galpy natural units
    R = (
        R_kpc
        /
        RO_MODEL
    )


    # Initially cold radial velocities
    vR = np.zeros(
        n
    )


    # Circular velocity in the background MW potential
    vT = vcirc(

        MWPotential2014,

        R,

        use_physical=False,
    )


    # 2D galpy Orbit initial conditions:
    #
    # [R, vR, vT, phi]

    return np.column_stack(

        [
            R,
            vR,
            vT,
            phi,
        ]
    )


# ============================================================
# 11. SIMULATE ONE PATTERN SPEED
# ============================================================

def simulate(
    pattern_speed,
    initial_conditions
):

    pot = make_spiral_potential(
        pattern_speed
    )

    orbit = Orbit(

        initial_conditions,

        ro=RO_MODEL,

        vo=VO_MODEL,
    )


    times = (

        np.linspace(

            0.0,

            T_END,

            N_OUTPUT,
        )

        *
        u.Gyr
    )


    # --------------------------------------------------------
    # FAST C INTEGRATOR
    # --------------------------------------------------------
    #
    # leapfrog_c sacrifices some accuracy compared with
    # higher-order methods, but is extremely fast.

    orbit.integrate(

        times,

        pot,

        method="leapfrog_c",

        dt=DT * u.Gyr,

        progressbar=False,
    )


    # Final state
    R = np.asarray(

        orbit.R(
            times[-1],
            use_physical=True
        )
    )


    phi = np.asarray(

        orbit.phi(
            times[-1],
        )
    )


    vR = np.asarray(

        orbit.vR(
            times[-1],
            use_physical=True
        )
    )


    vT = np.asarray(

        orbit.vT(
            times[-1],
            use_physical=True
        )
    )


    return (
        R,
        phi,
        vR,
        vT,
    )


# ============================================================
# 12. RUN ALL FAST SIMULATIONS
# ============================================================

def run_simulations():

    print(
        "\n"
        "========================================\n"
        "Running fast galpy spiral simulations\n"
        "========================================\n"
    )

    initial_conditions = initial_cold_disc()
    

    results = {}


    for omega in PATTERN_SPEEDS:

        print(

            f"omega_p = {omega:5.1f} "
            "km/s/kpc"
        )


        results[omega] = simulate(

            omega,

            initial_conditions
        )


    return results, initial_conditions


# ============================================================
# 13. FIGURE 4 — SIMULATED STELLAR ORBITS
# ============================================================

def plot_simulated_orbits(initial_conditions):

    selected = [

        PATTERN_SPEEDS[0],

        PATTERN_SPEEDS[
            len(PATTERN_SPEEDS) // 2
        ],

        PATTERN_SPEEDS[-1],
    ]


    fig, axes = plt.subplots(

        1,

        3,

        figsize=(14, 4.5)
    )


    # Using identical initial conditions for all three panels.

    initial = initial_conditions[:8]


    times = (

        np.linspace(

            0.0,

            T_END,

            N_OUTPUT,
        )

        *
        u.Gyr
    )


    for ax, omega in zip(
        axes,
        selected
    ):

        pot = make_spiral_potential(
            omega
        )


        orbit = Orbit(

            initial,

            ro=RO_MODEL,

            vo=VO_MODEL,
        )


        orbit.integrate(

            times,

            pot,

            method="leapfrog_c",

            dt=DT * u.Gyr,

            progressbar=False,
        )


        for j in range(
            8
        ):

            R = np.asarray(
                orbit.R(times, use_physical=True)
            )[j]

            phi = np.asarray(
                orbit.phi(times)
            )[j]


            x = (
                R
                *
                np.cos(phi)
            )

            y = (
                R
                *
                np.sin(phi)
            )


            ax.plot(

                x,

                y,

                linewidth=1.0,
            )


        # Solar radius marker
        ax.scatter(

            RO_MODEL,

            0,

            marker="*",

            s=100,
        )


        ax.set_xlim(-11.0, 11.0)
        ax.set_ylim(-11.0, 11.0)


        ax.set_aspect(
            "equal"
        )


        ax.set_xlabel(
            "x [kpc]"
        )

        ax.set_ylabel(
            "y [kpc]"
        )


        ax.set_title(

            rf"$\omega_p={omega:.0f}$ "
            r"km s$^{-1}$ kpc$^{-1}$"
        )


    fig.suptitle(

        "Stellar orbital response to different spiral pattern speeds"
    )


    fig.tight_layout()


    fig.savefig(

        "figure4_simulated_orbits.png",

        dpi=400,

        bbox_inches="tight",
    )


    plt.show()


# ============================================================
# 14. FIGURE 5 — PATTERN SPEED SCAN
# ============================================================

def plot_pattern_speed_scan(
    results,
    gaia,
):

    # --------------------------------------------------------
    # Gaia m=2 signal
    # --------------------------------------------------------

    gaia_amplitude, gaia_phase = (

        m2_fourier(

            gaia["R_kpc"].to_numpy(),

            gaia["phi_rad"].to_numpy(),

            gaia["vR_kms"].to_numpy(),
            phi_center=np.pi,
            phi_halfwidth=PHI_WEDGE_HALFWIDTH,
        )
    )


    # --------------------------------------------------------
    # Model response
    # --------------------------------------------------------

    model_amplitudes = []

    model_phases = []


    for omega in PATTERN_SPEEDS:

        R, phi, vR, vT = (
            results[omega]
        )


        amplitude, phase = (

            m2_fourier(

                R,

                phi,

                vR,
                phi_center=np.pi,
                phi_halfwidth=PHI_WEDGE_HALFWIDTH,
            )
        )


        model_amplitudes.append(
            amplitude
        )

        model_phases.append(
            phase
        )


    model_amplitudes = np.asarray(
        model_amplitudes
    )


    model_phases = np.asarray(
        model_phases
    )

    valid = np.isfinite(model_amplitudes)
    if not np.any(valid):
        raise RuntimeError(
            "No finite model m=2 amplitudes were measured. Increase N_SIM_STARS "
            "or widen PHI_WEDGE_HALFWIDTH."
        )

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    fig, ax = plt.subplots(

        figsize=(8, 5)
    )


    ax.plot(

        PATTERN_SPEEDS[valid],

        model_amplitudes[valid],

        marker="o",

        label="Exploratory galpy model",
    )


    if np.isfinite(
        gaia_amplitude
    ):

        ax.axhline(

            gaia_amplitude,

            linestyle="--",

            label=(

                rf"Gaia DR3: "
                rf"$A_2={gaia_amplitude:.2f}$ "
                r"km s$^{-1}$"
            ),
        )


    ax.set_xlabel(

        r"Spiral pattern speed "
        r"$\omega_p$ "
        r"[km s$^{-1}$ kpc$^{-1}$]"
    )


    ax.set_ylabel(

        r"$m=2$ radial-velocity "
        r"amplitude [km s$^{-1}$]"
    )


    ax.set_title(

        "Exploratory pattern-speed response scan"
    )


    ax.legend()


    fig.tight_layout()


    fig.savefig(

        "figure5_pattern_speed_scan.png",

        dpi=400,

        bbox_inches="tight",
    )


    plt.show()


    print()
    print(
        "Gaia DR3 m=2 radial-velocity amplitude:"
    )
    print(
        gaia_amplitude
    )

    if np.isfinite(gaia_amplitude):
        best_idx = np.nanargmin(np.abs(model_amplitudes - gaia_amplitude))
        print(
            f"Closest exploratory model response: omega_p = {PATTERN_SPEEDS[best_idx]:.1f} "
            f"km/s/kpc (A2 = {model_amplitudes[best_idx]:.2f} km/s)"
        )

    print()
    print(
        "IMPORTANT:"
    )

    print(
        "The pattern-speed comparison is exploratory."
    )

    print(
        "It is NOT a statistically complete inference "
        "of the Milky Way spiral pattern speed."
    )

def load_gaia():

    if GAIA_FILE.exists():

        print(
            f"Loading cached Gaia catalogue: {GAIA_FILE}"
        )

        return pd.read_csv(GAIA_FILE)

    try:
        return query_gaia()

    except Exception as exc:

        raise RuntimeError(
            "\nGaia DR3 could not be queried at the moment.\n"
            "This is usually an archive/TAP service issue rather "
            "than a problem with the simulation code.\n"
            "Retry later, or place a previously downloaded Gaia "
            "CSV at:\n"
            f"    {GAIA_FILE}\n\n"
            f"Underlying error:\n{exc}"
        )
# ============================================================
# 15. MAIN PROGRAM
# ============================================================

def main():

    print(
        "\n"
        "====================================================\n"
        " Gaia DR3 + galpy Milky Way Spiral-Arm Project\n"
        "====================================================\n"
    )


    # --------------------------------------------------------
    # Gaia
    # --------------------------------------------------------

    raw_gaia = load_gaia()


    gaia = transform_gaia(
        raw_gaia
    )


    print()


    # --------------------------------------------------------
    # GAIA FIGURES
    # --------------------------------------------------------

    print(
        "\nGenerating Figure 1..."
    )

    plot_gaia_radial_velocity_map(
        gaia
    )


    print(
        "Generating Figure 2..."
    )

    plot_gaia_velocity_distribution(
        gaia
    )


    print(
        "Generating Figure 3..."
    )

    plot_gaia_azimuthal_profile(
        gaia
    )


    # --------------------------------------------------------
    # GALPY SIMULATIONS
    # --------------------------------------------------------

    results, initial_conditions = run_simulations()


    # --------------------------------------------------------
    # SIMULATION FIGURE
    # --------------------------------------------------------

    print(
        "\nGenerating simulation Figure 4..."
    )

    plot_simulated_orbits(initial_conditions)


    # --------------------------------------------------------
    # PATTERN-SPEED SCAN
    # --------------------------------------------------------

    print(
        "Generating pattern-speed Figure 5..."
    )

    plot_pattern_speed_scan(

        results,

        gaia,
    )


    print(
        "\n"
        "====================================================\n"
        "DONE\n"
        "====================================================\n"
    )


if __name__ == "__main__":

    main()