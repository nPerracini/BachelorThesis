import os

from functions import bin_spectra
from functions import bin_axis
from functions import read_json
from functions import load_corrections
from functions import load_wavelengths
from functions import  load_data
from functions import  deskew_data

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# create spectrum of
def create_spectrum(folder, file_name, windows ,plot = True, show_plot = False, save_plot = False, corrected_chip = True, single = False):
    science_file = os.path.join("Data/" + folder, file_name)

    # load corrections
    bias, darkrate, gain, skews, rotation, beam_start, beam_end = load_corrections()

    # load data
    corr_data, header, obs_time, times = load_data(science_file, bias, darkrate, gain)

    # rotate data
    deskewed_data = deskew_data(corr_data, rotation, beam_start, beam_end, skews)

    # create wavelength array
    spectrum_raw = deskewed_data.sum(axis=0)[15:-15]
    spectrum = spectrum_raw

    # correct chip for wavelength bias
    if(corrected_chip):
        factors = read_json("detector.json")
        for i in range(len(spectrum)):
            spectrum[i] = spectrum[i] * factors[i]

    #load wavelenghts
    wavelengths = load_wavelengths(spectrum)

    # bin 2
    spectrum = bin_spectra(spectrum)
    wavelengths = bin_axis(spectrum, wavelengths)
    times = bin_axis(spectrum, times)

    # create alpha and constant windows, 935 and 1000 are pixel values around alpha and constant window
    H_alpha = []
    H_alpha_wavelength = []
    constant_zone_wavelength = []
    constant_zone = []
    for window in windows:
        H_alpha.append(spectrum[935-int(window/2) : 935+int(window/2)])
        H_alpha_wavelength.append(wavelengths[935-int(window/2) : 935+int(window/2)])

        constant_zone_wavelength.append(wavelengths[1000-int(window/2) : 1000+int(window/2)])
        constant_zone.append(spectrum[1000-int(window/2) : 1000+int(window/2)])


    # plot all spectra
    if(plot and not single):
        fig = plt.figure(num=2)
        plt.plot(wavelengths, spectrum, drawstyle="steps-mid", linewidth=0.5)
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("Intensity")
        plt.title("Spectrum of " + file_name)
        if(save_plot):
            plt.savefig(os.path.abspath(os.getcwd()) + "/Plots/" + folder + "/" + file_name + '.png')
            plt.clf()
        if(show_plot):
            plt.show(block=True)

    # plot single spectrum
    if (plot and single):
        # plot whole spectrum
        fig = plt.figure(num=2)
        plt.plot(wavelengths, spectrum, drawstyle="steps-mid", linewidth=0.5)
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("Intensity")
        plt.title("Spectrum of " + file_name)
        if (save_plot):
            plt.savefig(os.path.abspath(os.getcwd()) + "/Plots/" + folder + "/" + file_name + '.png')
            plt.clf()
        if (show_plot):
            plt.show(block=True)

        # plot focused spectrum
        fig = plt.figure(num=2)
        plt.plot(H_alpha_wavelength[2], H_alpha[2], linewidth=0.5)
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("Intensity")
        plt.title("Focused spectrum of " + file_name)
        if (save_plot):
            plt.savefig(os.path.abspath(os.getcwd()) + "/Plots/Focused_Spectrum.png")
            plt.clf()
        if (show_plot):
            plt.show(block=True)

    return H_alpha, constant_zone, H_alpha_wavelength, constant_zone_wavelength, times, spectrum_raw


# prepare data to plot alpha and constant windows using different windowsizes
def flux_graph(alpha_list, constant_list):
    alpha_flux = []
    constant_flux = []
    normal_flux = []
    for i in range(int(len(alpha_list))):
        for j in range(len(windows)):
            a = sum(alpha_list[i][j])
            c = sum(constant_list[i][j])
            alpha_flux.append(a)
            constant_flux.append(c)
            normal_flux.append(a / c)

    flux_plot(times, [alpha_flux[0::3], constant_flux[0::3]], [f"Alpha Line", f"Reference Window"], f"Summed Intensity Window Size {windows[0]}")
    flux_plot(times, [alpha_flux[1::3], constant_flux[1::3]], [f"Alpha Line", f"Reference Window"], f"Summed Intensity Window Size {windows[1]}")
    flux_plot(times, [alpha_flux[2::3], constant_flux[2::3]], [f"Alpha Line", f"Reference Window"], f"Summed Intensity Window Size {windows[2]}")
    flux_plot(times, [normal_flux[0::3], normal_flux[1::3], normal_flux[2::3]], ["1.70 nm Window", "4.25 nm Window", "8.50 nm Window"], "Normalized Intensity of H-Alpha line")

# plot all spectra above one another
def stacked_spectra_plot(spectra):
    wavelengths = load_wavelengths(spectra[0])
    for i in range(len(spectra)):
        spectra[i] = bin_spectra(spectra[i])
    wavelengths = bin_axis(spectra[0], wavelengths)

    plt.figure(figsize=(8, 10))
    for i in range(len(spectra)):
        spectra[i] += i*(np.nanmedian(spectra[0])/2)
        plt.plot(wavelengths, spectra[i], linewidth=0.3)
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Stacked Intensities")
    plt.title("Stacked Spectra")
    # add ticks on right side if needed
    #plt.yticks([spectra[22][len(spectra[22]) - 1], spectra[32][len(spectra[32]) - 1], spectra[40][len(spectra[40]) - 1]],["23", "33", "41"])
    plt.yticks([])
    plt.gca().yaxis.tick_right()
    plt.show(block=True)

# plot alpha around possible flares
def flare_plot(spectra, wavelengths, starting_index):
    plt.figure(figsize=(8, 5))
    for i in range(len(spectra)):
        plt.plot(wavelengths[i], spectra[i], label=f"Measurement {starting_index + i + 1}")
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Intensity")
    plt.title("H Alpha Line around possible Flare")
    plt.legend(loc=3)
    plt.show(block=True)

# prepare data for calculating and plotting, adjust plot, save_plot, show_plot and corrected_chip as needed
def prepare(folder, spectrum, windows):
    alpha, constant, alpha_wavelength, constant_wavelength, times, spectrum = create_spectrum(folder, spectrum, windows, plot=False,
                                                                             save_plot=False, show_plot=False,
                                                                             corrected_chip=True, single=single_spectrum)
    return alpha, alpha_wavelength, constant, constant_wavelength, times, spectrum


# plot intensity
def flux_plot(times, fluxes, legend, title):
    fig = plt.figure(num=2)
    for flux in fluxes:
        plt.plot(times, flux, linewidth=1)
    plt.gcf().autofmt_xdate()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.legend(legend, loc=3)
    plt.xlabel("Time [UTC]")
    plt.ylabel("Intensity")
    plt.title(title)
    #plt.ylim(ymin=0)
    plt.show(block=True)

if __name__ == "__main__":
    # check if only one spectrum is investigated
    single_spectrum = False
    # insert data folder
    folder = "AD_Leonis_new"
    spectra = []
    alpha_list = []
    alpha_wavelengths = []
    constant_list = []
    constant_wavelengths = []
    times = []
    # lenghts of windows
    windows = [10, 25, 50]
    counter = 0

    # prepare sinlge data
    if single_spectrum:
        alpha, alpha_wavelength, constant, constant_wavelength, t, spectrum = prepare(folder, "AD_Leo_100ym_300s_20260430_194408-9.fit", windows)

    # prepare full data set
    else:
        # sort by date
        files = sorted(
            os.listdir("Data/" + folder),
            key=lambda f: os.path.getctime(os.path.join("Data/" + folder, f))
        )

        # cut files if only looking at specific spectra, include first index, excludes last index
        # AD_Leonis_new: 9:14 24:29, 30:35, 38:43,   20:25
        # AD_Leonis: 0:4, 7:12, 10:15 ,13:18
        # AD_Leonis_100ym: 21:26, 25:30
        start = 10
        end = 15
        files = files[start:end]
        # read alpha and constant zones, and times
        for spectrum in files:
            alpha, alpha_wavelength, constant, constant_wavelength, t, spec = prepare(folder, spectrum, windows)
            spectra.append(spec)
            alpha_list.append(alpha)
            alpha_wavelengths.append(alpha_wavelength)
            constant_list.append(constant)
            constant_wavelengths.append(constant_wavelength)
            times.append(t)

        a = []
        aw = []
        c = []
        cw = []
        # prepare data to plot around flare, index indicates window size
        for i in range(len(alpha_list)):
            a.append(alpha_list[i][0])
            aw.append(alpha_wavelengths[i][0])
            c.append(constant_list[i][0])
            cw.append(constant_wavelengths[i][0])

        # plot alpha around flare
        flare_plot(a, aw, start)

    # plot stacked spectra
    stacked_spectra_plot(spectra)
    # plot alpha and constant
    flux_graph(alpha_list, constant_list)