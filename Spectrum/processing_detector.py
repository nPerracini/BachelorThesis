import os

import scipy.signal

from functions import linear
from functions import save_json
from functions import load_corrections
from functions import deskew_data_master

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import json


# load the data of an object
def get_image_data(object):
    files = []
    deskewed_data = []

    # load files
    for spectra in os.listdir("Data/" + object):
        files.append(os.path.join("Data", object + "/" + spectra))

    # load all corrections
    bias, darkrate, gain, skews, rotation, beam_start, beam_end = load_corrections()

    # apply all corrections
    for file in files:
        deskewed_data = deskew_data_master(file, bias, darkrate, gain, rotation, beam_start, beam_end, skews, deskewed_data)

    # prepare data
    deskewed_data = deskewed_data / len(files)
    spectrum = deskewed_data.sum(axis=0)[15:-15]

    # create wavelength array
    with open("calibration.json", "r") as f:
        calibration = json.load(f)

    lamda_0 = calibration["lamda_0"]
    lamda_1 = calibration["lamda_1"]

    wavelengths = []
    for i in range(len(spectrum)):
        wavelengths.append(linear(i, lamda_0, lamda_1))

    return wavelengths, spectrum

# create factors to calibrate the detector wavelength bias
if __name__ == "__main__":
    measured_wavelength, measured_intensity = get_image_data("Vega")

    # file of standard
    fits_file = "Data/Standard/alpha_lyr_mod_005.fits"

    # load standard
    with fits.open(fits_file) as hdul:
        raw_data = hdul[1].data

    # read wavelengths and intensity of standard
    vega_wavelength = []
    vega_intensity = []
    for item in raw_data:
        text_data = str(item).split(",")
        angstrom = float(text_data[0].strip("("))
        if(angstrom > 4500 and angstrom < 7500):
            vega_wavelength.append(angstrom/10)
            vega_intensity.append(float(text_data[1]))

    indices = []

    # interpolate standard data
    interpolated_flux = np.interp(measured_wavelength, vega_wavelength, vega_intensity)

    # calculate the correction factors
    factors = interpolated_flux/measured_intensity

    # apply smoothing filters
    filtered_factors_50 = scipy.signal.savgol_filter(factors, 50, 2)
    filtered_factors_150 = scipy.signal.savgol_filter(factors, 150, 2)
    filtered_factors_300 = scipy.signal.savgol_filter(factors, 300, 2)

    # adjust measured intensity
    original_intensity = []
    for i in range(len(measured_intensity)):
        original_intensity.append(measured_intensity[i]*factors[0])

    # correct measured intensity
    corrected_intensity = []
    for i in range(len(measured_intensity)):
        corrected_intensity.append(measured_intensity[i]*factors[i])

    # save to file
    save_json(filtered_factors_300.tolist(), "detector.json")

    # plot spectra
    plt.figure(figsize=(8, 5))
    plt.plot(vega_wavelength, vega_intensity, drawstyle="steps-mid", linewidth=1, color="r", label="Standard")
    plt.plot(measured_wavelength, original_intensity, drawstyle="steps-mid", linewidth=1, color="g", label="Measurement")
    plt.plot(measured_wavelength, corrected_intensity, drawstyle="steps-mid", linewidth=0.5, color="b", label="Correced Measurement")
    plt.legend()
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Intensity [erg/s cm^2 Å]")
    plt.title("Corrected Vega Spectrum")
    plt.show()

    # plot factors
    plt.figure(figsize=(8, 5))
    plt.plot(measured_wavelength, factors, drawstyle="steps-mid", linewidth=1, label="Raw Factors")
    plt.plot(measured_wavelength, filtered_factors_50, drawstyle="steps-mid", linewidth=1, label="Windowsize 50")
    plt.plot(measured_wavelength, filtered_factors_150, drawstyle="steps-mid", linewidth=1, label="Windowsize 150")
    plt.plot(measured_wavelength, filtered_factors_300, drawstyle="steps-mid", linewidth=1, label="Windowsize 300")
    plt.legend()
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Factor")
    plt.title("Detector Bias")
    plt.show()
