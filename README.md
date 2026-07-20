# Spectroscopic Measurement of Flare Stars

This code is part of a Bachelor Thesis by N. Perracini, 2026.
For access to observations, contact the N. Perracini.

## LightCurve
This sections explains the files in the LightCurve folder.  
The following structure is expected:
+ LightCurve
  - brightness_curve.py
  - flux_measurement.py
  - master.py
  - AD_LEO
  - Bias
  - Dark
  - master_bias.fits
  - master_dark_bias_corrected.fits
  - master_flat.tif

### brightness_curve.py
This code calculates the flux of AD Leonis and a reference star for a folder containing .fits observation files.  
It includes a correction for bias, dark and flats which require the following files:
      "master_bias.fits", "master_dark_bias_corrected.fits", "master_flat.tif"
+ The code provides 2 plots:
  - X-Axis: Time, Y-Axis: Intensity; This plot contains both AD Leonis and referenc star data and shows a brightness                 curve.
  - X-Axis: Time, Y-Axis: Intensity; This plot contain the normalized intensity of AD Leonis and shows a brightness                  curve.
This file needs hardcoded values for where the stars are in the image, in order to find the centroid.  
Adjust the name of the input folder if needed.  

### flux_measurement.py
This code illustrates the flux calculation using a single exposure of AD Leonis.  
It includes a correction for bias, dark and flats which require the following files:
      "master_bias.fits", "master_dark_bias_corrected.fits", "master_flat.tif"  
+ The code provides 1 plot:
  - X-Axis: pixel, Y-Axis: pixel; This plot contains the exposure. Marked are both stars, the red square is searched                 for the centroid, the inner annulus contains the area over which intensity is summed over, the band between                  annulus 2 and 3 contains the area over which the background radiation is measured.
This file needs hardcoded values for where the stars are in the image, in order to find the centroid.  
Adjust the name of the input file if needed.  

### master.py
This code creates master dark and bias frames. Use folders "Bias" and "Dark".  


## Spectrum
This sections explains the files in the Spectrum folder.  
The following structure is expected:
+ Spectrum
  - brightness_curve_lisa.py
  - functions.py
  - processing_detector.py
  - processing_full_data.py
  - processing_master_spectra.py
  - processing_neon.py
  - calibration.json
  - detector.json
  - Data
    - AD_Leonis
    - Bias
    - Dark
    - Capella
    - Flat
    - Jupiter
    - Neon
    - Standard
    - Vega
  - Plots
    - AD_Leonis
  - Processing_Data
    - Average_Bias_20260318.fits
    - Average_Dark_20260318.fits
    - Average_Dark_Rate_20260318.fits
    - Average_SkyFlat_20260318.fits
    - beam_edges.txt
    - beam_rotation.txt
    - Processed_Gaintable_20260318.fits
    - Slit_Skews_20260318.fits

 ### brightness_curve_lisa.py
 This code summs up the intensity of spectra inside a folder. The spectra are corrected for bias, dark and flats, as well as  adjusted for hot pixels and rotation.  
 + The code provides 1 plot:
   - X-Axis: Time, Y-Axis: Intensity; This plot shows the changing intensity over time and shows a brightness curve.
 Adjust the name of the input folder if needed.  

 ### functions.py
 This code includes functions used in the other files.  

 ### processing_detector.py
 This code corrects spectra of a standard star for bias, dark and flats, as well as adjusts them for hot pixels and rotation. The spectrum of the standard star is loaded and cut down to the necessary wavelengths. Factors to correct the detector bias are calculated and smoothed over. The factors are saved to "detector.json".  
 + The code provides 1 plot:
   - X-Axis: Wavelength, Y-Axis: Factor; This plot shows the factor needed at each wavelength.
 Adjust the name of the input folder if needed.  
 Adjust the name of the standard star if needed.  

 ### processing_full_data.py
 This code creates the spectrum for each exposure inside a folder. The spectra are corrected for bias, dark and flats, as well as adjusted for hot pixels and rotation. The spectra are binned over, combining 2 datapoints. Windows of different sizes are put over the H-Alpha line and a reference line.  
 + The code allows 6 plots:
   - X-Axis: Wavelength, Y-Axis: Intensity; This plot shows the spectrum of a measurement. (Inside create_spectrum                  function)
   - X-Axis: Wavelength, Y-Axis: Intensity; This plot shows a zoomed spectrum around the H-Alpha line of a measurement.             (Inside create_spectrum function)
   - X-Axis: Time, Y-Axis: Intensity; This plot shows the intensity summed over the H-Alpha line and a reference window             for different window sizes. (flux_plot function)
   - X-Axis: Time, Y-Axis: Intensity; This plot shows the normalized intensity over the H-Alpha line for different                  window sizes. (flux_plot function)
   - X-Axis: Wavelength, Y-Axis: Time/Intensity; This plot shows all spectra plotted above one another.                             (stacked_spectra_plot function)
   - X-Axis: Wavelength, Y-Axis: Intensity; This plot shows the H-Alpha intensity for a collection of measurements. The             index of the measurements need to be hardcoded. (flare_plot function)
Adjust name of the input folder if needed.  
Change single_spectrum variable if only one spectrum should be analyzed.  
Adjust window sizes and measurement indices if needed.  

### processing_master_spectra.py
This code creates the spectrum for each exposure inside a folder. The spectra are corrected for bias, dark and flats, as well as adjusted for hot pixels and rotation. The spectra are added up to create a master spectra.  
+ The code provides 2 plots:
  - X-Axis: pixel, Y-Axis: pixel; This plot shows a corrected and added up frame of the exposures.
  - X-Axis: Wavelength, Y-Axis: Intensity; This plot shows the master spectrum of the object.
Adjust name of input folder if needed.  

### processing_neon.py
This code creates the spectrum for each neon exposure inside a folder. The spectra are corrected for bias, dark and flats, as well as adjusted for hot pixels and rotation. The spectra are added up to create a master spectra. A peak detection finds the emission lines which are compared to literature data. A linear fit assigns every pixel value a wavelength. The wavelengths are saved to "calibration.json".  
+ The code provides 2 plots:
  - X-Axis: pixel, Y-Axis: pixel; This plot shows a corrected and added up frame of the exposures.
  - X-Axis: Wavelength, Y-Axis: Intensity; This plot shows the spectrum of the neon exposure. The detected peaks are               marked with x.
Adjust name of input folder if needed.  
Literature values need to be hardcoded.  
Detected peaks not in literature need to be hardcode excluded.  

### calibration.json
Contains the wavelength calibration from processing_neon.py  

### detector.json
Contains the factors for detector wavelength bias from processing_detector.py  

## Issues
+ Processing in LightCurve/brightness_curve.py needs to be adjusted. Either provide small images or optimize the process to decrease time for processing the images.\
+ Detection of bright pixels in LightCurve/brightness_curve.py and LightCurve/flux_measurement.py can be optimized to not need hardcoding. The presence of other stars and rotation during the night make this difficult.\
+ Empty folders might lead to problems. Error handling needs expansion.\
+ For processing in Spectrum/ the correction files were provided by Dr. Sellers, contact N. Perracini for access.\
