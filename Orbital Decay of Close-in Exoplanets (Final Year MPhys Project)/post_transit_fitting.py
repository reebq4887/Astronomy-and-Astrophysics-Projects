import pypickle
import corner
import matplotlib.pyplot as plt
import numpy as np
import glob

#path = "Processed Kepler Data/pkl files/long cadence"

def load_emcee(pkl_path):
  with open(pkl_path, "rb") as f:
    obj = pypickle.load(pkl_path)
    if isinstance(obj, dict):
      if 'posterior_samples' in obj:
        return obj['posterior_samples']
      return obj
  raise ValueError(f"Unrecognised pickle structure in {pkl_path}: {type(obj)}")

instrument = 'TESS'

# PARAMS = ['P_p1', 't0_p1', 'b_p1','p_p1','rho',f'q1_{instrument}',f'q2_{instrument}', f'mflux_{instrument}', f'sigma_w_{instrument}']
PARAMS = ['P_p1', 't0_p1','b_p1','p_p1','rho',f'q1_{instrument}',f'q2_{instrument}']

pkl_files = glob.glob('C:\\Users\\DELL\\Documents\\Astrophysics 4th Year Notes\\PX452 - Astrophysics Final Year Project\\Kepler-1658b\\posterior distribution\\emcee_combined.pkl', recursive=True)

rng = np.random.default_rng(123)
K = 1000

combined = {p: [] for p in PARAMS}

for pkl in pkl_files:
  samples = load_emcee(pkl)

  if not all(p in samples for p in PARAMS):
    continue

  for p in PARAMS:
    # combined[p].append(subsample(samples[p], K, rng))
    combined[p].append(samples[p])

for p in PARAMS:
  combined[p] = np.concatenate(combined[p])

X = np.vstack([combined[p] for p in PARAMS]).T

print("X shape:", X.shape)          # (Nsamples, Nparams)
print("Nsamples:", X.shape[0])

fig = corner.corner(
    X,
    labels=[r"$P$", r"$t_0$", r"$b$", r"$R_p/R_\star$", r"$\rho_\star$", r"$q_1$", r"$q_2$"],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_fmt=".3g",
    title_kwargs={"fontsize": 20},
    max_n_ticks=0,
    hist_kwargs={"linewidth":1.0},
)

for ax in fig.get_axes():

    # increase label font
    ax.xaxis.label.set_fontsize(16)
    ax.yaxis.label.set_fontsize(16)

    # move labels closer to plots
    ax.xaxis.set_label_coords(0.5, -0.05)
    ax.yaxis.set_label_coords(-0.05, 0.5)

    ax.set_title(ax.get_title(), pad=4)

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

plt.show()