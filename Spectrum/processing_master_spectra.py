import os

from functions import linear
from functions import read_json
from functions import load_corrections
from functions import deskew_data_master

import matplotlib.pyplot as plt


# create master spectrum
def process_images(object, detector=True):
    files = []
    deskewed_data = []

    # load files
    for spectra in os.listdir("Data/" + object):
        files.append(os.path.join("Data", object + "/" + spectra))

    # load corrections
    bias, darkrate, gain, skews, rotation, beam_start, beam_end = load_corrections()

    # correct and rotate data
    for file in files:
        deskewed_data = deskew_data_master(file, bias, darkrate, gain, rotation, beam_start, beam_end, skews, deskewed_data)

    # prepare data
    deskewed_data = deskewed_data / len(files)
    spectrum = deskewed_data.sum(axis=0)[15:-15]

    # correct detector wavelength bias
    if detector:
        factors = read_json("detector.json")
        for i in range(len(spectrum)):
            spectrum[i] = spectrum[i] * factors[i]

    # create wavelength array
    calibration = read_json("calibration.json")

    lamda_0 = calibration["lamda_0"]
    lamda_1 = calibration["lamda_1"]

    wavelengths = []
    for i in range(len(spectrum)):
        wavelengths.append(linear(i, lamda_0, lamda_1))

    # plot corrected and processed exposure
    fig = plt.figure(num=1)
    plt.imshow(deskewed_data, origin="lower", cmap="gray")
    plt.title("Corrected " + object + " frame")
    plt.xlabel("x [px]")
    plt.ylabel("y [x]")
    plt.show(block=False)

    # plot spectrum
    fig = plt.figure(num=2)
    plt.plot(wavelengths, spectrum, drawstyle="steps-mid", linewidth=0.5)
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Intensity")
    plt.title("Spectrum of " + object)
    plt.ylim(ymin = 0)
    plt.show(block=True)


# create master spectra of provided object
if __name__ == "__main__":
    process_images("Jupiter", detector=True)
