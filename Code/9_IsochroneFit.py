from IsochroneFitFunctions import estimate, plot_bma_HR
from SedFitFunctions import mean_density_from_MR
from inlist import *
import pandas as pd
import os
import numpy as np
import matplotlib

matplotlib.use('Agg')  # Non-interactive backend; avoid opening GUI windows

# MIST bands available for fitting
mist_available_bands = ['J', 'H', 'K', 'G', 'BP', 'RP', 'W1', 'W2', 'W3', 'TESS', 'Kepler']

# Mapping: CSV band name -> MIST column name
mist_band_to_column = {
    'G': 'gaiag',
    'BP': 'gaiabp',
    'RP': 'gaiarp',
    'J': 'twomassj',
    'H': 'twomassh',
    'K': 'twomassk',
    'W1': 'wisew1',
    'W2': 'wisew2',
    'W3': 'wisew3',
    'TESS': 'tess',
    'Kepler': 'kepler'
}
# Auto-built map: lowercase MIST column name -> canonical band label
band_name_map = {v.lower(): k for k, v in mist_band_to_column.items()}

# Load photometry table and SED fit outputs
photo = pd.read_csv(file_input_photometry, dtype={"starID": str})
sed = pd.read_csv(combined_data_sedfit_path + "Output_SED_Fits_Final.csv", dtype={"starID": str})
bands_used_all = pd.read_csv(combined_data_sedfit_path + "Output_SED_Fits_bandsUsed.csv", dtype={"starID": str})

# Align photometry rows with SED table by starID
photo = pd.merge(sed[["starID"]], photo, on="starID").reset_index(drop=True)
if not (photo.starID.equals(sed.starID)):
    print("The tables of input photometry and blackbody fits are not consistent, exiting.")
    exit()

photo = photo[np.append("starID", photo.columns[photo.columns != "starID"].sort_values())].copy()

# Fit all stars
all_results = []
for index, row in sed.iterrows():
    star_id = row["starID"]
    photo_row = photo[photo["starID"] == star_id]
    band_row = bands_used_all[bands_used_all["starID"] == star_id]

    if photo_row.empty or band_row.empty:
        print(f"Missing data for starID {star_id}, skipping.")
        continue

    # Bands listed as used; keep only those supported by MIST
    raw_band_list = band_row["band"].values[0].split('|')
    usable_bands = []
    for b in raw_band_list:
        std_b = band_name_map.get(b.lower())
        if std_b in mist_available_bands:
            usable_bands.append(std_b)

    # Build photometric parameters
    params = {
        "Teff": (row["modelteff"], row["modeltefferr"]),
        "logg": (row["modellogg"], row["modelloggerr"]),
        "feh": (row["feh"], row["feherr"]),
        "distance": (row["d"], row["derr"])
    }
    photo_dict = {}

    for b in mist_available_bands:
        if b not in usable_bands:
            continue
        col_name = mist_band_to_column.get(b)
        if col_name is None:
            continue
        mag_col = col_name
        err_col = col_name + "_e"
        if mag_col in photo_row.columns and err_col in photo_row.columns:
            mag_val = photo_row[mag_col].values[0]
            err_val = photo_row[err_col].values[0]
            if pd.notnull(mag_val) and pd.notnull(err_val):
                if err_val == 0.0:
                    err_val = abs(mag_val) * 0.005
                params[b] = (mag_val, err_val)
                photo_dict[b] = (mag_val, err_val)
    
    # Replace zero/invalid uncertainties for all tuple parameters
    for key in params:
        if isinstance(params[key], tuple):
            val, err = params[key]
            if err == 0.0 or np.isnan(err):  # Avoid divide-by-zero / NaN in likelihood
                params[key] = (val, max(abs(val) * 0.005, 1e-3))  # Floor: 0.5% of |value| or 1e-3

    # Save bands and photometry used for this star's fit
    band_photodf = pd.DataFrame({
        "starID": [star_id] * len(photo_dict),
        "band": list(photo_dict.keys()),
        "mag": [photo_dict[b][0] for b in photo_dict],
        "mag_err": [photo_dict[b][1] for b in photo_dict]
    })
    band_output_file = os.path.join(single_data_isochronefit_path, f"{star_id}_BandsUsed.csv")
    band_photodf.to_csv(band_output_file, index=False)

    # Start `isochrone fitting`
    print("-" * 15, f"Start Isochrone Fitting For {star_id}", "-" * 15)
    try:
        #med_logg, logg_error = estimate(usable_bands, params, logg=True)
        age_samples, mass_samples, eep_samples, logL_samples, radius_samples = estimate(
            usable_bands, params, logg=False)
    except Exception as e:
        print(f"Error fitting star {star_id}: {e}")
        continue

    # MIST posterior summaries (logL: log10(L/Lsun); radius: R/Rsun)
    lumi_mean = float(np.mean(logL_samples))
    lumi_std = float(np.std(logL_samples))
    radius_mean = float(np.mean(radius_samples))
    radius_std = float(np.std(radius_samples))

    # Organize the fitting results (inputs Teff/logg/feh/d match SED columns used in params)
    results = {
        "starID": star_id,
        "Teff": params["Teff"][0],
        "Tefferr": params["Teff"][1],
        "logg": params["logg"][0],
        "loggerr": params["logg"][1],
        "feh": params["feh"][0],
        "feherr": params["feh"][1],
        "d": params["distance"][0],
        "derr": params["distance"][1],
        "lumi_MIST": lumi_mean,
        "lumierr_MIST": lumi_std,
        "radius_MIST": radius_mean,
        "radiuserr_MIST": radius_std,
        "Mass_MIST": float(np.mean(mass_samples)),
        "Masserr_MIST": float(np.std(mass_samples)),
        "EEP_MIST": float(np.mean(eep_samples)),
        "EEPerr_MIST": float(np.std(eep_samples)),
        "Age_MIST": float(np.mean(age_samples)),
        "Ageerr_MIST": float(np.std(age_samples)),
    }

    # Data output for plotting
    out = {
        'mist_samples': {
            'age': age_samples
        },
        'weighted_average': {
            'z': [params['feh'][0] for _ in range(len(age_samples))]
        }
    }

    # Save this star's isochrone fit results
    single_star_df = pd.DataFrame([results])
    star_output_file = os.path.join(single_data_isochronefit_path, f"{star_id}.csv")
    single_star_df.to_csv(star_output_file, index=False)

    plot_bma_HR(results, out, method='samples', nsamp=5, png=True, out_folder=fig_isochronefit_path)
    print(f"Saved results for {star_id} to {star_output_file}")
    print("-" * 25, "End Isochrone Fitting", "-" * 25)
    print("%" * 70)

    all_results.append(results)
    
    # print the results
    print(f"Used bands for {star_id}: {usable_bands}")
    print(f"Teff: {results['Teff']:.2f} ± {results['Tefferr']:.2f}, "
          f"logg: {results['logg']:.2f} ± {results['loggerr']:.2f}, "
          f"[Fe/H]: {results['feh']:.2f} ± {results['feherr']:.2f}, "
          f"d: {results['d']:.2f} ± {results['derr']:.2f}")
    print(f"log10(L/Lsun) [MIST]: {results['lumi_MIST']:.3f} ± {results['lumierr_MIST']:.3f}, "
          f"R/Rsun [MIST]: {results['radius_MIST']:.3f} ± {results['radiuserr_MIST']:.3f}, "
          f"Mass: {results['Mass_MIST']:.3f} ± {results['Masserr_MIST']:.3f}, "
          f"EEP: {results['EEP_MIST']:.3f} ± {results['EEPerr_MIST']:.3f}, "
          f"Age: {results['Age_MIST']:.3f} ± {results['Ageerr_MIST']:.3f}")

# Combined output for all successfully fitted stars
if all_results:
    combined_results_df = pd.DataFrame(all_results)
    combined_output_file = os.path.join(combined_data_isochronefit_path, "Output_Isochrone_Fits.csv")
    combined_results_df.to_csv(combined_output_file, index=False)
    print(f"Saved all results to {combined_output_file}")
else:
    print("No stars were successfully fitted, check for errors.")

# Combined SED + isochrone table (enable via inlist.Lcombine_sed_isochrone)
if Lcombine_sed_isochrone and all_results:
    iso_df = pd.DataFrame(all_results)
    sed_cols = [
        "starID", "Av", "Averr", "ScaleF", "ScaleFerr",
        "lumi", "lumierr", "radius", "radiuserr", "Fbol", "Fbolerr",
    ]
    for _opt in ("lumibol", "lumibolerr", "Redchi"):
        if _opt in sed.columns and _opt not in sed_cols:
            sed_cols.append(_opt)
    sed_cols = [c for c in sed_cols if c in sed.columns]
    sed_sub = sed[sed_cols].copy()
    comb = iso_df.merge(sed_sub, on="starID", how="left")
    rename_map = {
        "lumi": "lumi_SED",
        "lumierr": "lumierr_SED",
        "radius": "radius_SED",
        "radiuserr": "radiuserr_SED",
        "lumibol": "lumibol_SED",
        "lumibolerr": "lumibolerr_SED",
        "Redchi": "Redchi_SED",
    }
    rename_map = {k: v for k, v in rename_map.items() if k in comb.columns}
    comb = comb.rename(columns=rename_map)
    # rho_bar = 3M/(4 pi R^3) in g/cm^3; MIST uses MIST M+R; SED uses same Mass_MIST with SED radius
    _mist_pairs = [
        mean_density_from_MR(
            r["Mass_MIST"], r["radius_MIST"],
            r["Masserr_MIST"], r["radiuserr_MIST"],
        )
        for _, r in comb.iterrows()
    ]
    _sed_pairs = [
        mean_density_from_MR(
            r["Mass_MIST"], r["radius_SED"],
            r["Masserr_MIST"], r["radiuserr_SED"],
        )
        for _, r in comb.iterrows()
    ]
    comb["meandensity_MIST"], comb["meandensity_MIST_err"] = zip(*_mist_pairs)
    comb["meandensity_SED"], comb["meandensity_SED_err"] = zip(*_sed_pairs)
    column_order = [
        "starID", "Teff", "Tefferr", "logg", "loggerr", "feh", "feherr", "d", "derr",
        "Av", "Averr", "ScaleF", "ScaleFerr",
        "lumi_MIST", "lumierr_MIST", "radius_MIST", "radiuserr_MIST",
        "Mass_MIST", "Masserr_MIST", "EEP_MIST", "EEPerr_MIST", "Age_MIST", "Ageerr_MIST",
        "lumi_SED", "lumierr_SED", "radius_SED", "radiuserr_SED", "Fbol", "Fbolerr",
        "meandensity_MIST", "meandensity_MIST_err", "meandensity_SED", "meandensity_SED_err",
    ]
    extra = [c for c in ["lumibol", "lumibolerr", "Redchi_SED"] if c in comb.columns]
    out_cols = [c for c in column_order + extra if c in comb.columns]
    comb_out = comb[out_cols].copy()
    os.makedirs(os.path.dirname(file_output_combine) or ".", exist_ok=True)
    comb_out.to_csv(file_output_combine, index=False)
    print(f"Saved SED+isochrone combine table to {file_output_combine}")