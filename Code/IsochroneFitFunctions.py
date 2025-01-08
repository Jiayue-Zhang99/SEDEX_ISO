"""Estimate logg using MIST isochrones."""
"""isochrones packages can be installed from: https://isochrones.readthedocs.io/en/latest/install.html"""

import sys
import warnings
import pickle

import pandas as pd
import numpy as np
from isochrones import SingleStarModel, get_ichrone
from isochrones.mist import MIST_Isochrone
from isochrones.priors import FlatPrior, GaussianPrior
from numba.core.errors import (NumbaDeprecationWarning,
                               NumbaPendingDeprecationWarning)
import dynesty
from dynesty.utils import resample_equal

import traceback
from scipy.special import erf

from random import choice
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os

warnings.filterwarnings(
    'ignore', category=NumbaDeprecationWarning, append=True)
warnings.filterwarnings(
    'ignore', category=NumbaPendingDeprecationWarning, append=True)

#from error import (DynestyError, InputError)
class Error(Exception):
    """Base class for exceptions in this module."""

    def __repr__(self):
        """Error identification for logging."""
        return self.errorname

    def __str__(self):
        """Error identification for logging."""
        return self.errorname

    def __raise__(self):
        """Raise an exception and print the error message."""
        self.warn()
        sys.exit()

    def warn(self):
        """Print error message."""
        print('An exception was caught!', end=': ')
        print(self, end='\nError message: ')
        print(self.message)

    def log(self, out):
        """Log the error."""
        log_f = open(out, 'a')
        log_f.write(self.message)
        log_f.close()

    pass

class DynestyError(Error):
    """Exception raised when dynesty crashes."""

    def __init__(self, out, mod, e):
        self.errorname = 'DynestyError'
        self.message = 'ERROR OCCURRED DURING DYNESTY RUN.\n'
        self.message += f'WHILE FITTING MODEL {mod}\n'
        self.message += f'DUMPING `sampler.results` TO {out}\n'
        self.message += f'ERROR READS:\n{traceback.format_exc()}'

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes
    ----------
    message -- explanation of the error

    """

    def __init__(self, message):
        self.errorname = 'InputError'
        self.message = message

# from utils import credibility_interval
def credibility_interval(post, alpha=1.):
    """Calculate bayesian credibility interval.

    Parameters:
    -----------
    post : array_like
        The posterior sample over which to calculate the bayesian credibility
        interval.
    alpha : float, optional
        Confidence level.
    Returns:
    --------
    med : float
        Median of the posterior.
    low : float
        Lower part of the credibility interval.
    up : float
        Upper part of the credibility interval.

    """
    z = erf(alpha / np.sqrt(2))

    lower_percentile = 100 * (1 - z) / 2
    upper_percentile = 100 * (1 + z) / 2
    low, med, up = np.percentile(
        post, [lower_percentile, 50, upper_percentile]
    )
    return med, low, up

############### Functions for Isochrone Fitting ##############

def get_isochrone(logage, feh):
    """Retrieve isochrone for given age and feh."""
    mist = MIST_Isochrone()
    iso = mist.isochrone(logage, feh)
    return iso


def estimate(bands, params, logg=True, out_folder='.'):
    """Estimate logg using MIST isochrones."""
    mist = get_ichrone('mist', bands=bands)
    model = SingleStarModel(mist, **params)
    if 'distance' in params.keys():
        dist, dist_e = params['distance']
    elif 'parallax' in params.keys():
        dist = 1 / (params['parallax'][0] * 0.001)
        dist_e = dist * params['parallax'][1] / params['parallax'][0]
    else:
        msg = 'No parallax or distance found.'
        msg += 'Aborting age and mass calculation.'
        InputError(msg).warn()
        return np.zeros(10), np.zeros(10)
    if 'feh' in params.keys():
        fe, fe_e = params['feh']
        if fe + fe_e >= 0.5:
            model._priors['feh'] = FlatPrior([-0.5, 0.5])
        else:
            model._priors['feh'] = GaussianPrior(fe, fe_e)
    if 'mass' in params.keys():
        m, m_e = params['mass']
        model._priors['mass'] = GaussianPrior(m, m_e)
    if 'AV' in params.keys():
        av, av_e = params['AV']
        model._priors['AV'] = GaussianPrior(av, av_e)
    model._priors['distance'] = GaussianPrior(dist, dist_e)
    sampler = dynesty.NestedSampler(
        loglike, prior_transform, model.n_params + len(bands),
        nlive=500, bound='multi', sample='rwalk',
        logl_args=([model, params, bands]),
        ptform_args=([model])
    )
    try:
        sampler.run_nested(dlogz=0.01)
    except ValueError as e:
        dump_out = f'{out_folder}/isochrone_DUMP.pkl'
        pickle.dump(sampler.results, open(dump_out, 'wb'))
        DynestyError(dump_out, 'isochrone', e).__raise__()
    results = sampler.results
    samples = resample_equal(
        results.samples, np.exp(results.logwt - results.logz[-1])
    )
    ###########################################################################
    # Written by Dan Foreman-mackey
    # https://github.com/dfm/gaia-isochrones
    df = model._samples = pd.DataFrame(
        dict(
            zip(
                list(model.param_names),
                samples.T,
            )
        )
    )
    model._derived_samples = model.ic(
        *[df[c].values for c in model.param_names])
    model._derived_samples["parallax"] = 1000.0 / df["distance"]
    model._derived_samples["distance"] = df["distance"]
    model._derived_samples["AV"] = df["AV"]
    ###########################################################################
    if logg:
        samples = model._derived_samples['logg']
        med, lo, up = credibility_interval(samples, 5)
        med_e = max([med - lo, up - med])
        return med, med_e
    else:
        age_samples = 10 ** (model._derived_samples['age'] - 9)
        mass_samples = model._derived_samples['mass']
        eep_samples = model._derived_samples['eep']
        return age_samples, mass_samples, eep_samples


# Written by Dan Foreman-mackey
# https://github.com/dfm/gaia-isochrones

# These functions wrap isochrones so that they can be used with dynesty:
def prior_transform(u, mod):
    cube = np.copy(u)
    mod.mnest_prior(cube[: mod.n_params], None, None)
    cube[mod.n_params:] = -10 + 20 * cube[mod.n_params:]
    return cube


def loglike(theta, mod, params, jitter_vars):
    ind0 = mod.n_params
    lp0 = 0.0
    for i, k in enumerate(jitter_vars):
        err = np.sqrt(params[k][1] ** 2 + np.exp(theta[ind0 + i]))
        lp0 -= 2 * np.log(err)  # This is to fix a bug in isochrones
        mod.kwargs[k] = (params[k][0], err)
    lp = lp0 + mod.lnpost(theta[: mod.n_params])
    if np.isfinite(lp):
        return np.clip(lp, -1e10, np.inf)
    return -1e10


def plot_bma_HR(results, out, method, nsamp, hr_figsize=(10, 8), hr_cmap='plasma',
                hr_color='blue', hr_marker='o', fontsize=16,
                tick_labelsize=14, png=False, pdf=False, out_folder='.'):
    """
    Plot HR diagram for the star.
    """
    print(f'Plotting HR diagram for Star ID: {results["starID"]}')
    
    # Read stellar parameters
    starID = results["starID"]
    age = results["Age"]
    feh = results["feh"]
    teff = results["Teff"]
    lum = results["lumi"]
    
    # Compute error in log space to avoid horizontal/vertical error bars
    teff_err = np.abs(np.log10(teff + results["Tefferr"]) - np.log10(teff))
    lum_err = np.abs(np.log10(lum + results["lumierr"]) - np.log10(lum))
    
    # Convert Teff and Luminosity to log scale
    log_teff = np.log10(teff)
    log_lum = np.log10(lum)

    # Handle upper limit for FeH
    if feh > 0.5:
        feh = 0.5

    # Retrieve sampled isochrones from data
    ages = out['mist_samples']['age']
    m = method if method == 'samples' else 'average'
    fehs = out[f'weighted_average']['z']

    # Create the figure
    fig, ax = plt.subplots(figsize=hr_figsize, dpi=192)
    
    # === Step 1: Plot multiple sampled isochrones (gray) in the background ===
    for _ in range(nsamp):
        a = np.log10(choice(ages)) + 9  
        z = choice(fehs)  
        if z > 0.5:
            z = 0.5
        iso = get_isochrone(a, z)
        ax.plot(iso['logTeff'].values, iso['logL'].values, color='gray', alpha=0.5, linewidth=1.5, zorder=1)

    # === Step 2: Plot best-fitting isochrone on top ===
    iso_bf = get_isochrone(np.log10(age) + 9, feh)
    logteff = iso_bf['logTeff'].values
    loglum = iso_bf['logL'].values
    mass = iso_bf['mass'].values

    # Use LineCollection to plot the best-fitting isochrone
    points = np.array([logteff, loglum]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(mass.min(), mass.max())
    lc = LineCollection(segments, cmap=hr_cmap, norm=norm, linewidths=3, zorder=2)  # Ensure it's above the gray lines

    lc.set_array(mass)
    ax.add_collection(lc)
    cbar = fig.colorbar(lc, ax=ax, pad=0.01)
    cbar.set_label(r'$M_\odot$', rotation=270, fontsize=fontsize, labelpad=20)

    # === Step 3: Plot the observed star (error bars and data point) ===
    ax.errorbar(log_teff, log_lum, xerr=[[teff_err], [teff_err]],
                yerr=[[lum_err], [lum_err]], color=hr_color, capsize=3, fmt='o', zorder=3)
    ax.scatter(log_teff, log_lum, s=200, color=hr_color, edgecolors='k', marker=hr_marker, zorder=4, alpha=0.7)

    # Invert x-axis for HR diagram
    ax.invert_xaxis()
    ax.set_xlabel(r'$\log(T_{\rm eff})$', fontsize=fontsize)
    ax.set_ylabel(r'$\log L$', fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=tick_labelsize)
    ax.set_title(f'Star ID: {starID}', fontsize=fontsize+2)

    # Ensure the save path exists and save the figure
    os.makedirs(out_folder, exist_ok=True)
    save_path_png = os.path.join(out_folder, f"{starID}.png")
    save_path_pdf = os.path.join(out_folder, f"{starID}.pdf")

    if png:
        plt.savefig(save_path_png, bbox_inches='tight')
    if pdf:
        plt.savefig(save_path_pdf, bbox_inches='tight')

    plt.show()
    print(f"HR Diagram saved to {save_path_png}")