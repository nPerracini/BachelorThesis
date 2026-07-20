import os

from functions import linear
from functions import save_json
from functions import read_json
from functions import load_corrections
from functions import deskew_data_master

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy.optimize as optimize


# calculate wavelength calibration
def wavelength_calibration(peaks):
    # numbers 4, 9, 23, 24 are peaks not in literature
    clean_peaks = []
    for i in range(len(peaks)):
        if(i != 3 and i != 8 and i != 22 and i != 23):
            clean_peaks.append(peaks[i])

    # set literature peak wavelengths
    literature_peaks  = [585.249, 588.190, 594.483, 603.000, 607.434,
                         609.616, 614.306, 621.728, 626.650, 630.479,
                         633.443, 638.299, 640.225, 650.653, 653.288,
                         659.895, 667.828, 671.704, 692.947, 703.241]

    # create linear fit
    params, cov = optimize.curve_fit(linear, clean_peaks, literature_peaks)
    lamda_0, lamda_1 = params

    # print calibration
    errors = np.sqrt(np.diag(cov))
    print(f"Linear fit: lamda(x) = [{lamda_0} +- {errors[0]}] + [{lamda_1} +- {errors[1]} * x] ")

    # save to file
    calibration = {"lamda_0": float(lamda_0), "lamda_1": float(lamda_1),
                   "error_0": float(errors[0]), "error_1": float(errors[1])}
    save_json(calibration, "calibration.json")


# calculate wavelength calibration
if __name__ == "__main__":
    files = []
    deskewed_data = []

    # load exposures
    for spectra in os.listdir("Data/Neon"):
        files.append(os.path.join("Data", "Neon/" + spectra))

    # load corrections
    bias, darkrate, gain, skews, rotation, beam_start, beam_end = load_corrections()

    # apply corrections and rotations
    for file in files:
        deskewed_data = deskew_data_master(file, bias, darkrate, gain, rotation, beam_start, beam_end, skews, deskewed_data)

    # prepare data
    deskewed_data = deskewed_data / len(files)
    spectrum = deskewed_data.sum(axis=0)[15:-15]

    # load detector bias
    factors = read_json("detector.json")
    for i in range(len(spectrum)):
        spectrum[i] = spectrum[i] * factors[i]

    # find peaks and their width
    peaks, _ = signal.find_peaks(spectrum, prominence=0.02*np.max(spectrum), width=(1,8))
    widths, width_heights, left_ips, right_ips = signal.peak_widths(spectrum, peaks, rel_height=0.5)

    # plot corrected neon exposure
    fig = plt.figure(num=1)
    plt.imshow(deskewed_data, origin="lower", cmap="gray")
    plt.title("Corrected neon frame")
    plt.xlabel("x [px]")
    plt.ylabel("y [x]")
    plt.show(block=False)

    # plot neon spectrum with marked peaks
    fig = plt.figure(num=2)
    plt.plot(spectrum, drawstyle="steps-mid", linewidth=0.5)
    for i in range (len(peaks)):
        if(i != 3 and i != 8 and i != 22 and i != 23):
            plt.plot(peaks[i], spectrum[peaks[i]], marker="x", color="red", linestyle="None")
    plt.xlabel("x [px]")
    plt.ylabel("Intensity")
    plt.title("Spectrum of internal neon lamp")
    plt.show(block=True)

    # run calibration
    wavelength_calibration(peaks)
