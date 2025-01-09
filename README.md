# SEDEX_ISO: An Extended Version of SEDEX with Isochrone Fitting Capability

**SEDEX_ISO** is an enhanced version of SEDEX, incorporating the **isochrone fitting** functionalities from the [astroARIADNE](https://github.com/jvines/astroARIADNE.git) codebase. This integration utilizes the **ISOCHRONE package**, allowing SEDEX_ISO to provide a more comprehensive analysis pipeline for stellar parameters.

## Key Features of SEDEX_ISO

SEDEX_ISO extends the capabilities of SEDEX by adding an automated workflow that sequentially performs:
1. **Blackbody Radiation Fitting** – Initial estimation of stellar parameters.
2. **SED Fitting** – Utilizing the MARCS and BOSZ spectral libraries for a more precise determination of stellar properties.
3. **Isochrone Fitting** – Leveraging the **ISOCHRONE package** to derive stellar ages and evolutionary states.

Compared with the original SEDEX, SEDEX_ISO **adds isochrone-based age and mass estimations**, making it a powerful tool for studying stellar evolution.

## Input and Data Handling

SEDEX_ISO follows the core input structure of SEDEX but expands upon it:

- **Input:** A list of target stars, specified in `UserInputData.csv`, containing:
  - Gaia DR3 `source_id`
  - Stellar atmospheric parameters:  
    $T_{\rm{eff}}$, $\sigma_{T_{\rm{eff}}}$, $\log g$, $\sigma_{\log g}$, $\textrm{[Fe/H]}$, $\sigma_{\textrm{[Fe/H]}}$
  - Extinction parameters: $A_V$, $\sigma_{A_V}$ (optional)

- **Photometry:** SEDEX_ISO automatically retrieves broadband photometry from multiple surveys, including:
  - Gaia (`G`, `BP`, `RP`)
  - 2MASS (`J`, `H`, `K_S`)
  - WISE (`W1`, `W2`, `W3`, `W4`)
  - Hipparcos, Tycho2, SDSS, Pan-STARRS1, APASS, SkyMapper, and more.

- **Isochrone Models:** SEDEX_ISO reads the **stellar parameters derived from SED Fitting** (such as $T_{\rm{eff}}$, $\sigma_{T_{\rm{eff}}}$, $\log g$, $\sigma_{\log g}$, $\textrm{[Fe/H]}$, $\sigma_{\textrm{[Fe/H]}}$, Gaia photometry (`G`, `BP`, `RP`)) and uses them as inputs for **isochrone fitting**. The fitting process employs **isochrone grids** from the **ISOCHRONE package**, allowing SEDEX_ISO to determine the **stellar age, mass, and evolutionary stage** with improved accuracy.


## Output

SEDEX_ISO generates the following results:

- **SED Fitting Parameters** (same as SEDEX)
  - Effective temperature
  - Surface gravity
  - Extinction
- **Isochrone-Derived Properties:**
  - Stellar age
  - Mass
  - Evolutionary stage classification
- **Fitting Figures:**
  - Blackbody fitting
  - SED fitting
  - Isochrone fitting

## How to Use SEDEX_ISO

The usage of SEDEX_ISO closely follows SEDEX with additional isochrone fitting steps:

1. **Set up a Gaia account** and input credentials in `inlist.py`.
2. **Download spectral libraries**: [MARCS](https://marcs.astro.uu.se/) (standard composition) and [BOSZ](https://archive.stsci.edu/hlsp/bosz/search.php) (solar [C/M], resolution $R=20000$).
3. **Prepare the target list** in `UserInputData.csv`.
4. **Configure `inlist.py`**, ensuring that the **isochrone fitting step is enabled**.
5. **Run `submit.sh`** to process the targets through:
   - Blackbody fitting
   - SED fitting
   - Isochrone fitting
6. **Retrieve results** from:
   - `../Data/Output_Fits/"sample"/SEDFits/Output_SED_Fits_Final.csv` (SED results)
   - `../Data/Output_Fits/"sample"/IsochorneFits/Output_Isochrone_Fits.csv` (new isochrone-derived parameters)
   - Figures in `../Figures/"sample"/`

By incorporating **isochrone fitting**, SEDEX_ISO provides a more **comprehensive** approach to stellar characterization, making it particularly useful for studies involving stellar ages, mass distributions, and evolutionary states.

---