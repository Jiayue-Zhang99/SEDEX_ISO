from IsochroneFitFunctions import estimate, plot_bma_HR
from inlist import *
import pandas as pd
import os
import numpy as np
import matplotlib

matplotlib.use('TkAgg')

# Read SED and photometry data
photo = pd.read_csv(file_input_photometry, dtype={"starID": str})  
sed = pd.read_csv(combined_data_sedfit_path + "Output_SED_Fits_Final.csv", dtype={"starID": str})

# Ensure that `photo` and `sed` only contain matching `starID`
photo = pd.merge(sed[["starID"]], photo, on="starID").reset_index(drop=True)

if not (photo.starID.equals(sed.starID)): 
    print("The tables of input photometry and blackbody fits are not consistent, exiting.")
    exit()

# Standardize the format of `photo`
subs = np.where(photo.columns != "starID")[0]
photo = photo[np.append("starID", photo.columns[subs].sort_values())].copy()

# Store all fitting results in a DataFrame
all_results = []

# Iterate through each star in `sed`
for index, row in sed.iterrows():
    star_id = row["starID"]

    # Find the corresponding data for `starID` in `photo`
    photo_row = photo[photo["starID"] == star_id]
    if photo_row.empty:
        print(f"Warning: starID {star_id} not found in photometry data, skipping.")
        continue

    # Generate `params`
    params = {
        "Teff": (row["modelteff"], row["modeltefferr"]),
        "logg": (row["modellogg"], row["modelloggerr"]),
        "feh": (row["feh"], row["feherr"]),  
        "distance": (row["d"], row["derr"]),
        "G": (photo_row["gaiag"].values[0], photo_row["gaiag_e"].values[0]),
        "BP": (photo_row["gaiabp"].values[0], photo_row["gaiabp_e"].values[0]),
        "RP": (photo_row["gaiarp"].values[0], photo_row["gaiarp_e"].values[0]),
    }
    
    # Ensure a minimum error value of 0.005 times the parameter value if the error is zero
    for key in ["Teff", "logg", "feh", "distance", "G", "BP", "RP"]:
        if params[key][1] == 0.0:
            params[key] = (params[key][0], params[key][0] * 0.005)

    # Start `isochrone fitting`
    print("-" * 15, f"Start Isochrone Fitting For {star_id}", "-" * 15)
    try:
        med_logg, logg_error = estimate(bands, params, logg=True)
        age_samples, mass_samples, eep_samples = estimate(bands, params, logg=False)
    except Exception as e:
        print(f"Error fitting star {star_id}: {e}")
        continue

    # Compute log errors to avoid horizontal lines in plots
    teff_log = np.log10(params["Teff"][0])
    teff_err = np.abs(np.log10(params["Teff"][0] + params["Teff"][1]) - teff_log)
    
    lum_mean = mass_samples.mean()
    lum_err = np.abs(np.log10(lum_mean + mass_samples.std()) - np.log10(lum_mean))

    # Organize the fitting results
    results = {
        "starID": star_id,
        "Teff": params["Teff"][0],
        "Tefferr": params["Teff"][1],
        "logg": med_logg,
        "loggerr": logg_error,
        "feh": params["feh"][0],
        "feherr": params["feh"][1],
        "d": params["distance"][0],
        "derr": params["distance"][1],
        "lumi": lum_mean,
        "lumierr": lum_err,
        "Mass": mass_samples.mean(),
        "Masserr": mass_samples.std(),
        "EEP": eep_samples.mean(),
        "EEPerr": eep_samples.std(),
        "Age": age_samples.mean(),
        "Ageerr": age_samples.std(),
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

    # Save the single star CSV file
    single_star_df = pd.DataFrame([results])
    star_output_file = os.path.join(single_data_isochronefit_path, f"{star_id}.csv")
    single_star_df.to_csv(star_output_file, index=False)

    # Generate the HR diagram and save
    plot_bma_HR(results, out, method='samples', nsamp=5, png=True, out_folder=fig_isochronefit_path)

    # Display the saved CSV path
    print(f"Saved results for {star_id} to {star_output_file}")
    print("-" * 25, "End Isochrone Fitting", "-" * 25)
    print("%" * 70)

    # Append results to the DataFrame
    all_results.append(results)

# Combine all results into a single CSV file and save as `Output_Isochrone_Fits.csv`
if all_results:
    combined_results_df = pd.DataFrame(all_results)
    combined_output_file = os.path.join(combined_data_isochronefit_path, "Output_Isochrone_Fits.csv")
    combined_results_df.to_csv(combined_output_file, index=False)
    print(f"Saved all results to {combined_output_file}")

else:
    print("No stars were successfully fitted, check for errors.")