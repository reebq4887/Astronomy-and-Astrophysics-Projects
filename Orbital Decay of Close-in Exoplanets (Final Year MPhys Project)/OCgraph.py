import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import glob
import os
import pandas as pd
from scipy.stats import chi2

# -------------------------
# Paths
# -------------------------
path_tess = "tess\*.dat"
path_kepler = "kepler\*.dat"

path = r"C:\Users\DELL\Documents\Astrophysics 4th Year Notes\PX452 - Astrophysics Final Year Project\Kepler-1658b\Processed Kepler Data\posteriors\global"

path_tess = glob.glob(
    os.path.join(
        path, path_tess
    )
)
path_kepler = glob.glob(
    os.path.join(
        path, path_kepler
    )
)

# -------------------------
# Constants
# -------------------------
#Constants
kepler_int = 2454833.0
tess_int = 2457000.0
G = 6.67430e-11
M_star = 1.45 * 1.989e30
s_Mstar = 0.06 * 1.989e30
M_planet = 5.88 * 1.898e27
s_Mplanet = 0.47 * 1.898e27
R_star = 2.9* 6.957e8
s_Rstar = 0.12 * 6.957e8
P_ref = 3.8493733
tref_bjd = 2455005.92415

a_over_R_tess = 4.1138766916
# a_over_R_kep = 4.7497460635
a_over_R_kep = 4.04

a_over_R_tess_plus = 0.1126432795
# a_over_R_kep_plus = 0.0932231456
a_over_R_kep_plus = 0.18
a_over_R_tess_minus = 0.1185534272
# a_over_R_kep_minus = 0.1098543582
a_over_R_kep_minus = 0.17

a_over_R_correct = 4.04
a_over_R_correct_plus = 0.18
a_over_R_correct_minus = 0.17

# x_ref = 2454833.0
# t_plot = t_bjd - x_ref

# -------------------------
# Helper functions
# -------------------------
def read_posteriors(file):
    return pd.read_csv(
        file,
        delim_whitespace=True,
        comment="#",
        names=["param", "median", "upper", "lower"],
    )

def get_t0_and_bounds(df):
    row = df[df["param"] == "t0_p1"].iloc[0]
    t0_median = row["median"]
    t0_plus = row["upper"]
    t0_minus = row["lower"]
    t0_upper = t0_median + t0_plus
    t0_lower = t0_median - t0_minus
    return t0_median, t0_plus, t0_minus, t0_upper, t0_lower

def computeoc(t0_bjd, P_ref, tref_bjd):
    E = np.rint((t0_bjd - tref_bjd) / P_ref).astype(int)
    tcalc = tref_bjd + E * P_ref
    oc_days = t0_bjd - tcalc
    return E, oc_days * 24 * 60  # minutes

def asym_summary(samples):
    q16, q50, q84 = np.percentile(samples, [16, 50, 84])
    return q50, q84 - q50, q50 - q16

def sample_asym(med, sig_plus, sig_minus, size=1, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    z = rng.normal(size=size)
    return np.where(z >= 0, med + z * sig_plus, med + z * sig_minus)


def weighted_chi2(E, t, sig, order=1):
    E = np.asarray(E)
    t = np.asarray(t)
    sig = np.asarray(sig)

    if order == 1:
        X = np.vstack([np.ones_like(E), E]).T
    elif order == 2:
        X = np.vstack([np.ones_like(E), E, E**2]).T
    else:
        raise ValueError("order must be 1 or 2")

    w = 1.0 / (sig**2)
    Aw = X * np.sqrt(w[:, None])
    yw = t * np.sqrt(w)

    coeff = np.linalg.lstsq(Aw, yw, rcond=None)[0]
    t_model = X @ coeff

    chi2_val = np.sum(((t - t_model) / sig) ** 2)
    dof = len(t) - len(coeff)
    return coeff, chi2_val, dof, t_model
# -------------------------
# Load TESS posteriors
# -------------------------
tess_t0bjd = []
tess_sig_plus = []
tess_sig_minus = []
tess_upper = []
tess_lower = []

kepler_t0bjd = []
kepler_sig_plus = []
kepler_sig_minus = []
kepler_upper = []
kepler_lower = []

for pf in path_tess:
    df = read_posteriors(pf)
    t0_median, t0_plus, t0_minus, t0_upper, t0_lower = get_t0_and_bounds(df)

    tess_t0bjd.append(t0_median + tess_int)
    tess_sig_plus.append(t0_plus)
    tess_sig_minus.append(t0_minus)
    tess_upper.append(t0_upper + tess_int)
    tess_lower.append(t0_lower + tess_int)

for pf in path_kepler:
    df = read_posteriors(pf)
    t0_median, t0_plus, t0_minus, t0_upper, t0_lower = get_t0_and_bounds(df)

    kepler_t0bjd.append(t0_median + kepler_int)
    kepler_sig_plus.append(t0_plus)
    kepler_sig_minus.append(t0_minus)
    kepler_upper.append(t0_upper + kepler_int)
    kepler_lower.append(t0_lower + kepler_int)

tess_t0bjd = np.array(tess_t0bjd)
sort = np.argsort(tess_t0bjd)
tess_t0bjd = tess_t0bjd[sort]
tess_sig_plus = np.array(tess_sig_plus)[sort]
tess_sig_minus = np.array(tess_sig_minus)[sort]
tess_upper = np.array(tess_upper)[sort]
tess_lower = np.array(tess_lower)[sort]

print(len(tess_t0bjd))

print("TESS t0 (BJD) =", tess_t0bjd)
print("TESS sigma plus (days) =", tess_sig_plus)
print("TESS sigma minus (days) =", tess_sig_minus)

kep_outlier = kepler_t0bjd[0]
kep_outlier_plus = kepler_sig_plus[0]
kep_outlier_minus = kepler_sig_minus[0]
kep_outlier_upper = kepler_upper[0]
kep_outlier_lower = kepler_lower[0]

tess_outlier = tess_t0bjd[1]
tess_outlier_plus = tess_sig_plus[1]
tess_outlier_minus = tess_sig_minus[1]
tess_outlier_upper = tess_upper[1]
tess_outlier_lower = tess_lower[1]

tess_t0bjd = np.delete(np.array(tess_t0bjd), 1)
tess_sig_plus = np.delete(np.array(tess_sig_plus), 1)
tess_sig_minus = np.delete(np.array(tess_sig_minus), 1)
tess_upper = np.delete(np.array(tess_upper), 1)
tess_lower = np.delete(np.array(tess_lower), 1)



kepler_t0bjd = np.array(kepler_t0bjd)[1:]
kepler_sig_plus = np.array(kepler_sig_plus)[1:]
kepler_sig_minus = np.array(kepler_sig_minus)[1:]
kepler_upper = np.array(kepler_upper)[1:]
kepler_lower = np.array(kepler_lower)[1:]

pal_wirc = np.array([2459097.8002, 2459790.6819])
pal_sig_plus = np.array([0.0015, 0.0015])
pal_sig_minus = np.array([0.0015, 0.0013])
pal_upper = np.array([2459097.8002+0.0015, 2459790.6819+0.0015])
pal_lower = np.array([2459097.8002-0.0015, 2459790.6819-0.0013])

# symmetric sigma only for weighting the least-squares fit
tess_sig = 0.5 * (tess_sig_plus + tess_sig_minus)
kep_sig = 0.5 * (kepler_sig_plus + kepler_sig_minus)
pal_sig = 0.5 * (pal_upper - pal_lower)

# -------------------------
# Compute O-C values
# -------------------------
E_kep, oc_kep_min = computeoc(kepler_t0bjd, P_ref, tref_bjd)
E_tess, oc_tess_min = computeoc(tess_t0bjd, P_ref, tref_bjd)
E_pal, oc_pal_min = computeoc(pal_wirc, P_ref, tref_bjd)

E_kep_outlier, oc_kep_outlier = computeoc(kep_outlier, P_ref, tref_bjd)
E_tess_outlier, oc_tess_outlier = computeoc(tess_outlier, P_ref, tref_bjd)

# Combined timings used in ephemeris fit
t_bjd = np.concatenate([tess_t0bjd, kepler_t0bjd, pal_wirc])
sig = np.concatenate([tess_sig, kep_sig, pal_sig])

sort = np.argsort(t_bjd)
t_bjd = t_bjd[sort]
sig = sig[sort]

E = np.rint((t_bjd - tref_bjd) / P_ref).astype(int)
w = 1.0 / np.maximum(sig, 1e-12) ** 2

# -------------------------
# Linear ephemeris fit
# -------------------------
A_lin = np.vstack([np.ones_like(E), E]).T
Aw_lin = A_lin * np.sqrt(w[:, None])
yw = t_bjd * np.sqrt(w)

coeff_lin = np.linalg.lstsq(Aw_lin, yw, rcond=None)[0]
a_lin, b_lin = coeff_lin

cov_lin = np.linalg.inv(A_lin.T @ (w[:, None] * A_lin))

t0_const = a_lin
P_const = b_lin

t0_const_err = np.sqrt(cov_lin[0, 0])
P_const_err = np.sqrt(cov_lin[1, 1])

print("\nConstant-period model parameters")
print(f"t0 (BJD_TDB) = {t0_const:.8f} ± {t0_const_err:.8f}")
print(f"P (days)     = {P_const:.10f} ± {P_const_err:.10f}")

# -------------------------
# Quadratic ephemeris fit
# -------------------------
A_quad = np.vstack([np.ones_like(E), E, E**2]).T
Aw_quad = A_quad * np.sqrt(w[:, None])

coeff_quad = np.linalg.lstsq(Aw_quad, yw, rcond=None)[0]
a, b, c = coeff_quad

cov_quad = np.linalg.inv(A_quad.T @ (w[:, None] * A_quad))

t0_decay = a
P_decay = b
dPdN = 2 * c   # IMPORTANT

# Uncertainties from covariance matrix
t0_decay_err = np.sqrt(cov_quad[0, 0])
P_decay_err = np.sqrt(cov_quad[1, 1])
c_err = np.sqrt(cov_quad[2, 2])
dPdN_err = 2 * c_err

print("\nOrbital decay model parameters")
print(f"t0 (BJD_TDB) = {t0_decay:.8f} ± {t0_decay_err:.8f}")
print(f"P (days)     = {P_decay:.10f} ± {P_decay_err:.10f}")
print(f"dP/dN (days/orbit) = {dPdN:.3e} ± {dPdN_err:.3e}")

Pdot = dPdN / P_decay   # days/day

# Convert to ms/yr
msday = 86400 * 1000
days_per_year = 365.25

Pdot_msyr = Pdot * msday * days_per_year
Pdot_msyr_err = (dPdN_err / P_decay) * msday * days_per_year

print(f"Pdot (ms/yr) = {Pdot_msyr:.3f} ± {Pdot_msyr_err:.3f}")



# ----------------------
# Posterior-sampled derived quantities
# -------------------------
N = 20000
rng = np.random.default_rng(42)
Nmc = 20000

draws_quad = rng.multivariate_normal([a, b, c], cov_quad, size=N)
a_s = draws_quad[:, 0]
b_s = draws_quad[:, 1]
c_s = draws_quad[:, 2]

Pdot_s = 2 * c_s / b_s

M_star_s = rng.normal(M_star, s_Mstar, N)
M_planet_s = rng.normal(M_planet, s_Mplanet, N)
R_star_s = rng.normal(R_star, s_Rstar, N)

# asymmetric sampling for a/R*
u = rng.normal(size=N)
a_over_R_kep_s = np.where(
    u >= 0,
    a_over_R_kep + u * a_over_R_kep_plus,
    a_over_R_kep + u * a_over_R_kep_minus
)

a_over_R_tess_s = np.where(
    u >= 0,
    a_over_R_tess + u * a_over_R_tess_plus,
    a_over_R_tess + u * a_over_R_tess_minus
)

a_over_R_correct_s= np.where(
    u >= 0,
    a_over_R_correct + u * a_over_R_correct_plus,
    a_over_R_correct + u * a_over_R_correct_minus
)

# Derived quantities

Pdot_mc = []
timescale_mc = []
for _ in range(Nmc):
    tess_draw = np.array([
        sample_asym(m, sp, sm, size=1, rng=rng)[0]
        for m, sp, sm in zip(tess_t0bjd, tess_sig_plus, tess_sig_minus)
    ])

    kepler_draw = np.array([
        sample_asym(m, sp, sm, size=1, rng=rng)[0]
        for m, sp, sm in zip(kepler_t0bjd, kepler_sig_plus, kepler_sig_minus)
    ])

    pal_draw = np.array([ sample_asym(m, sp, sm, size=1, rng=rng)[0] for m, sp, sm in zip(pal_wirc, pal_sig_plus, pal_sig_minus) ])

    t_draw = np.concatenate([tess_draw, kepler_draw, pal_draw])
    sig_draw = np.concatenate([
        0.5 * (tess_sig_plus + tess_sig_minus),   # fit weights
        0.5 * (kepler_sig_plus + kepler_sig_minus),
        0.5 * (pal_sig_plus + pal_sig_minus)
    ])

    sort = np.argsort(t_draw)
    t_draw = t_draw[sort]
    sig_draw = sig_draw[sort]

    E = np.rint((t_draw - tref_bjd) / P_ref).astype(int)
    w = 1.0 / np.maximum(sig_draw, 1e-12)**2

    A_quad = np.vstack([np.ones_like(E), E, E**2]).T
    Aw = A_quad * np.sqrt(w[:, None])
    yw = t_draw * np.sqrt(w)

    coeff = np.linalg.lstsq(Aw, yw, rcond=None)[0]
    a_mc, b_mc, c_mc = coeff

    Pdot_val = 2 * c_mc / b_mc
    tau_val = P_ref / np.abs(Pdot_val) / 365.25 / 1e6

    Pdot_mc.append(Pdot_val)
    timescale_mc.append(tau_val)

Pdot_mc = np.array(Pdot_mc)
timescale_mc = np.array(timescale_mc)

Pdot_med, Pdot_plus, Pdot_minus = asym_summary(Pdot_mc)
print(f"Pdot (days/day) = {Pdot_med:.3e} +{Pdot_plus:.3e} -{Pdot_minus:.3e}")

msday = 86400 * 1000.0
days_per_year = 365.25
Pdot_msyr_mc = Pdot_mc * msday * days_per_year

Pdot_msyr_med, Pdot_msyr_plus, Pdot_msyr_minus = asym_summary(Pdot_msyr_mc)
print(f"Pdot (ms/yr) mc = {Pdot_msyr_med:.3f} +{Pdot_msyr_plus:.3f} -{Pdot_msyr_minus:.3f}")

tau_med, tau_plus, tau_minus = asym_summary(timescale_mc)
print(f"timescale (Myr) mc = {tau_med:.3f} +{tau_plus:.3f} -{tau_minus:.3f}")



Pdot_med, Pdot_plus, Pdot_minus = asym_summary(Pdot_s)
# print(f"Pdot (days/day) = {Pdot_med:.3e} +{Pdot_plus:.3e} -{Pdot_minus:.3e}")

msday = 86400 * 1000.0
days_per_year = 365.25
Pdot_msyr_s = Pdot_s * msday * days_per_year
Pdot_msyr_med, Pdot_msyr_plus, Pdot_msyr_minus = asym_summary(Pdot_msyr_s)
print(f"Pdot (ms/yr) s = {Pdot_msyr_med:.3f} +{Pdot_msyr_plus:.3f} -{Pdot_msyr_minus:.3f}")

timescale_s = P_ref / np.abs(Pdot_s) / 365.25 / 1e6
tau_med, tau_plus, tau_minus = asym_summary(timescale_s)
print(f"timescale (Myr) s = {tau_med:.3f} +{tau_plus:.3f} -{tau_minus:.3f}")

# dist_s = a_over_R_s * R_star_s
# dist_AU_s = dist_s / 1.496e11
# dist_med, dist_plus, dist_minus = asym_summary(dist_AU_s)
# print(f"semimajor axis (AU) = {dist_med:.5f} +{dist_plus:.5f} -{dist_minus:.5f}")

# adot_s = (2.0 / 3.0) * (Pdot_s / P_ref) * dist_s / 86400.0
# adot_med, adot_plus, adot_minus = asym_summary(adot_s)
# print(f"adot (m/s) = {adot_med:.3e} +{adot_plus:.3e} -{adot_minus:.3e}")

# orbital_L_s = M_planet_s * np.sqrt(G * M_star_s * dist_s)
# L_med, L_plus, L_minus = asym_summary(orbital_L_s)
# print(f"Orbital Angular Momentum (kg m^2/s) = {L_med:.3e} +{L_plus:.3e} -{L_minus:.3e}")

# orbital_energy_s = -G * M_star_s * M_planet_s / (2 * dist_s)
# E_med, E_plus, E_minus = asym_summary(orbital_energy_s)
# print(f"Orbital Energy (J) = {E_med:.3e} +{E_plus:.3e} -{E_minus:.3e}")

# dldt_s = 0.5 * (adot_s / dist_s) * orbital_L_s
# dEdt_s = 0.5 * (adot_s / dist_s) * orbital_energy_s

# dL_med, dL_plus, dL_minus = asym_summary(dldt_s)
# dE_med, dE_plus, dE_minus = asym_summary(dEdt_s)

# print(f"dL/dt (kg m^2/s^2) = {dL_med:.3e} +{dL_plus:.3e} -{dL_minus:.3e}")
# print(f"dE/dt (W) = {dE_med:.3e} +{dE_plus:.3e} -{dE_minus:.3e}")

tidalq_kep_s = -27 * np.pi / (2 * Pdot_s) * (M_planet_s / M_star_s) * (1 / a_over_R_kep_s) ** 5
Q_med_kep, Q_plus_kep, Q_minus_kep = asym_summary(tidalq_kep_s)
print(f"Tidal Q Kepler s = {Q_med_kep:.3e} +{Q_plus_kep:.3e} -{Q_minus_kep:.3e}")

tidalq_tess_s = -27 * np.pi / (2 * Pdot_s) * (M_planet_s / M_star_s) * (1 / a_over_R_tess_s) ** 5
Q_med_tess, Q_plus_tess, Q_minus_tess = asym_summary(tidalq_tess_s)
print(f"Tidal Q TESS s = {Q_med_tess:.3e} +{Q_plus_tess:.3e} -{Q_minus_tess:.3e}")


tidalq_kep_s = -27 * np.pi / (2 * Pdot_mc) * (M_planet_s / M_star_s) * (1 / a_over_R_kep_s) ** 5
Q_med_kep, Q_plus_kep, Q_minus_kep = asym_summary(tidalq_kep_s)
print(f"Tidal Q Kepler mc= {Q_med_kep:.3e} +{Q_plus_kep:.3e} -{Q_minus_kep:.3e}")

tidalq_tess_s = -27 * np.pi / (2 * Pdot_mc) * (M_planet_s / M_star_s) * (1 / a_over_R_tess_s) ** 5
Q_med_tess, Q_plus_tess, Q_minus_tess = asym_summary(tidalq_tess_s)
print(f"Tidal Q TESS mc= {Q_med_tess:.3e} +{Q_plus_tess:.3e} -{Q_minus_tess:.3e}")


tidalq_correct = -27 * np.pi / (2 * Pdot_mc) * (M_planet_s / M_star_s) * (1 / a_over_R_correct_s) ** 5 * (2/3 - 1)
Q_med_correct, Q_plus_correct, Q_minus_correct = asym_summary(tidalq_correct)
print(f"Tidal Q CORRECT= {Q_med_correct:.3e} +{Q_plus_correct:.3e} -{Q_minus_correct:.3e}")

# -------------------------
# Smooth grid for O-C models
# -------------------------
tg = np.linspace(t_bjd.min(), t_bjd.max() + 400, 1200)
Eg = (tg - tref_bjd) / P_ref
trefgrid = tref_bjd + Eg * P_ref

# Best-fit model curves
tlingrid = a_lin + b_lin * Eg
tqgrid = a + b * Eg + c * Eg**2

oc_const = (tlingrid - trefgrid) * 24 * 60
oc_dec = (tqgrid - trefgrid) * 24 * 60

# Posterior-sampled uncertainty bands for the ephemerides
draws_lin_band = rng.multivariate_normal([a_lin, b_lin], cov_lin, size=5000)
models_lin = draws_lin_band[:, 0, None] + draws_lin_band[:, 1, None] * Eg[None, :]
oc_models_lin = (models_lin - trefgrid[None, :]) * 24 * 60
oc_const_lo = np.percentile(oc_models_lin, 16, axis=0)
oc_const_hi = np.percentile(oc_models_lin, 84, axis=0)

draws_quad_band = rng.multivariate_normal([a, b, c], cov_quad, size=5000)
models_quad = (
    draws_quad_band[:, 0, None]
    + draws_quad_band[:, 1, None] * Eg[None, :]
    + draws_quad_band[:, 2, None] * Eg[None, :]**2
)
oc_models_quad = (models_quad - trefgrid[None, :]) * 24 * 60
oc_dec_lo = np.percentile(oc_models_quad, 16, axis=0)
oc_dec_hi = np.percentile(oc_models_quad, 84, axis=0)

# -------------------------
# Plot styling
# -------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 14,
    "axes.labelsize": 22,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
})

fig, (ax, ax_res) = plt.subplots(
    2, 1, figsize=(12, 8.0), sharex=True,
    gridspec_kw={"height_ratios": [3.2, 1], "hspace": 0.05}
)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")
ax_res.set_facecolor("white")

# Shaded model bands
ax.fill_between(
    tg, oc_const_lo, oc_const_hi,
    color="cornflowerblue", alpha=0.55, linewidth=0, zorder=1
)
ax.fill_between(
    tg, oc_dec_lo, oc_dec_hi,
    color="peru", alpha=0.55, linewidth=0, zorder=1
)

# Main model curves
ax.plot(tg, oc_const, color="royalblue", lw=2.5,
        label="Linear ephemeris (constant period)", zorder=3)
ax.plot(tg, oc_dec, color="tomato", lw=2.5,
        label=r"Quadratic ephemeris (orbital decay, $\dot{P} < 0$)", zorder=3)

# Data points
ax.errorbar(
    kep_outlier,
    oc_kep_outlier,
    yerr=np.vstack([kep_outlier_minus * 24 * 60, kep_outlier_plus * 24 * 60]),
    fmt="o",
    ms=6.0,
    color="white",
    ecolor="tab:blue",
    elinewidth=1.5,
    capsize=3,
    markeredgecolor="tab:blue", markeredgewidth=0.4, zorder=4
)

ax.errorbar(
    tess_outlier,
    oc_tess_outlier,
    yerr=np.vstack([tess_outlier_minus * 24 * 60, tess_outlier_plus * 24 * 60]),
    fmt="o",
    ms=6.5,
    color="white",
    ecolor="orangered",
    elinewidth=1.5,
    capsize=3,
    markeredgecolor="orangered", markeredgewidth=0.4, zorder=4
)

ax.errorbar(
    tess_t0bjd,
    oc_tess_min,
    yerr=np.vstack([tess_sig_minus * 24 * 60, tess_sig_plus * 24 * 60]),
    fmt="o",
    ms=6.5,
    color="orangered",
    ecolor="orangered",
    elinewidth=1.5,
    capsize=3,
    markeredgecolor="black", markeredgewidth=0.4,
    label="TESS data", zorder=4
)
ax.errorbar(
    kepler_t0bjd,
    oc_kep_min,
    yerr=kep_sig * 24 * 60,
    fmt="o",
    ms=6.0,
    color="tab:blue",
    ecolor="tab:blue",
    elinewidth=1.5,
    capsize=3,
    markeredgecolor="black", markeredgewidth=0.4,
    label="Kepler data", zorder=4
)

ax.errorbar(
    pal_wirc,
    oc_pal_min,
    yerr=pal_sig * 24 * 60,
    fmt="o",
    ms=6.0,
    color="green",
    ecolor="green",
    elinewidth=1.5,
    capsize=3,
    markeredgecolor="black", markeredgewidth=0.4,
    label="Palomar/WIRC data (Vissapragada et al. 2022)", zorder=4
)

# Reference line
ax.axhline(0, color="gray", ls="--", lw=1.5, alpha=0.9)
ax_res.axhline(0, color="gray", ls="--", lw=1.2, alpha=0.9)

# Labels
ax.set_ylabel("Transit Timing Deviation (min)", fontsize=18)
ax_res.set_xlabel(r"BJD$_{\mathrm{TDB}}$ - 2457000 (days)",fontsize=18)
ax_res.set_ylabel("Residuals",fontsize=18)

# Limits and annotations on shifted x-axis
xmin_plot = (t_bjd.min()) - 80
xmax_plot = (t_bjd.max()) + 180
ax.set_xlim(xmin_plot, xmax_plot)
ax_res.set_xlim(xmin_plot, xmax_plot)
ax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
ax_res.xaxis.set_major_locator(ticker.MultipleLocator(1000))

# Residuals relative to quadratic model
oc_quad_tess = np.interp(tess_t0bjd, tg, oc_dec)
oc_quad_kep = np.interp(kepler_t0bjd, tg, oc_dec)
oc_quad_pal = np.interp(pal_wirc, tg, oc_dec)
res_tess = oc_tess_min - oc_quad_tess
res_kep = oc_kep_min - oc_quad_kep
res_pal = oc_pal_min - oc_quad_pal

oc_quad_kep_outlier = np.interp(kep_outlier, tg, oc_dec)
oc_quad_tess_outlier = np.interp(tess_outlier, tg, oc_dec)
res_kep_outlier = oc_kep_outlier - oc_quad_kep_outlier
res_tess_outlier = oc_tess_outlier - oc_quad_tess_outlier

ax_res.errorbar(
    kep_outlier, res_kep_outlier,
    yerr=np.vstack([kep_outlier_minus * 24 * 60, kep_outlier_plus * 24 * 60]),
    fmt="o", ms=5, color="white", ecolor="tab:blue",
    elinewidth=1.1, capsize=3, markeredgecolor="black", markeredgewidth=0.35,
    zorder=4
)

ax_res.errorbar(
    tess_outlier, res_tess_outlier,
    yerr=np.vstack([tess_outlier_minus * 24 * 60, tess_outlier_plus * 24 * 60]),
    fmt="o", ms=5, color="white", ecolor="orangered",
    elinewidth=1.1, capsize=3, markeredgecolor="black", markeredgewidth=0.35,
    zorder=4
)

ax_res.errorbar(
    tess_t0bjd, res_tess,
    yerr=np.vstack([tess_sig_minus * 24 * 60, tess_sig_plus * 24 * 60]),
    fmt="o", ms=5, color="orangered", ecolor="orangered",
    elinewidth=1.1, capsize=3, markeredgecolor="black", markeredgewidth=0.35,
    zorder=4
)
ax_res.errorbar(
    kepler_t0bjd, res_kep, yerr=kep_sig * 24 * 60,
    fmt="o", ms=5, color="tab:blue", ecolor="tab:blue",
    elinewidth=1.1, capsize=3, markeredgecolor="black", markeredgewidth=0.35,
    zorder=4
)
ax_res.errorbar(
    pal_wirc, res_pal, yerr=pal_sig * 24 * 60,
    fmt="o", ms=5, color="green", ecolor="green",
    elinewidth=1.1, capsize=3, markeredgecolor="black", markeredgewidth=0.35,
    zorder=4
)


ax.text(
    (t_bjd.max()),
    -2.8,
    "Chontos et al. (2019) Ephemerides",
    color="gray",
    fontsize=16,
    ha="right",
    alpha=0.9,
)

ax.legend(loc="lower left", frameon=True, fancybox=True, framealpha=0.9)
ax.grid(alpha=0.15)
ax_res.grid(alpha=0.15)

for spine in ax.spines.values():
    spine.set_linewidth(1.2)
for spine in ax_res.spines.values():
    spine.set_linewidth(1.2)

ax.ticklabel_format(style="plain", axis="x", useOffset=False)
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{x - tess_int:.0f}"))

plt.tight_layout()
plt.show()

# -------------------------
# Model comparison statistics
# -------------------------
theta_lin, chi2_lin, dof_lin, t_model_lin = weighted_chi2(E, t_bjd, sig, order=1)
theta_quad, chi2_quad, dof_quad, t_model_quad = weighted_chi2(E, t_bjd, sig, order=2)

print("Linear Fit: chi2 =", chi2_lin, ", dof =", dof_lin, ", chi2/dof =", chi2_lin / dof_lin)
print("Quadratic Fit: chi2 =", chi2_quad, ", dof =", dof_quad, ", chi2/dof =", chi2_quad / dof_quad)

delta_chi2 = chi2_lin - chi2_quad
p_value = 1 - chi2.cdf(delta_chi2, df=dof_lin - dof_quad)
print("Delta chi2 =", delta_chi2)
print("p-value for improvement from linear to quadratic fit:", p_value)

n = len(t_bjd)
klin, kquad = 2, 3
AIC_lin = chi2_lin + 2 * klin
AIC_quad = chi2_quad + 2 * kquad
BIC_lin = chi2_lin + klin * np.log(n)
BIC_quad = chi2_quad + kquad * np.log(n)

print("\nAIC:", AIC_lin, AIC_quad, "AIC difference:", AIC_lin - AIC_quad)
print("BIC:", BIC_lin, BIC_quad, "BIC difference:", BIC_lin - BIC_quad)