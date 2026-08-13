import numpy as np
import glob
import os
import pandas as pd

path_kep = glob.glob(os.path.join("Processed Kepler Data", "posteriors", "global","kepler", "*.dat"))
path_tess = glob.glob(os.path.join("Processed Kepler Data", "posteriors", "global","tess", "*.dat"))

kepler_int = 2454833.0
tess_int = 2457000.0

def read_posteriors(file):
    df = pd.read_csv(file, delim_whitespace=True, comment='#', 
                 names=['param', 'median', 'upper', 'lower'])
    return df

def get_param_and_sig(df, param):
    row = df[df['param'] == param].iloc[0]
    param_median = row['median']
    param_sig = 0.5 * (row['upper'] + row['lower'])
    return param_median, param_sig

def weighted_mean(param_median, param_sig):
    w = 1.0 / (param_sig**2)
    mu = np.sum(w * param_median) / np.sum(w)
    mu_err = np.sqrt(1.0 / np.sum(w))
    return mu, mu_err

instrument = 'TESS'
PARAMS = ['P_p1', 't0_p1','b_p1','p_p1','rho',f'q1_{instrument}',f'q2_{instrument}', f'mflux_{instrument}', f'sigma_w_{instrument}']


vals = []
sigs = []
for param in PARAMS:
    for pf in path_tess:
        df = read_posteriors(pf)
        param_median, param_sig = get_param_and_sig(df, param)
        vals.append(param_median)
        sigs.append(param_sig)
    vals = np.array(vals)
    sigs = np.array(sigs)
    mu, mu_err = weighted_mean(vals, sigs)
    with open('posterior_dist_tess.txt', 'a') as f:
        line = f'{param}'+' '+f'{mu}'+' '+f'{mu_err}\n'
        f.writelines(line)
    vals = []
    sigs = []


