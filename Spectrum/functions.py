import numpy as np
import astropy.io.fits as fits
import scipy.ndimage as ndimage
from astropy.time import Time
import json

"""This function was provided by Dr. Sellers"""
# remove hot pixels and cosmic rays
def image_despiker(image: np.ndarray, footprint: tuple=(5, 1), srange: tuple=(0.75, 1.25)) -> np.ndarray:
    """Simple median-filter based image despiker. Removes hot pixels and cosmic rays,
    backfills those pixels with median-filter values.

    Parameters
    ----------
    image : np.ndarray
        Image to despike
    footprint : tuple, optional
        Passed to scipy.ndimage.median_filter's size keyword
        Default (5, 1)
    srange : tuple, optional
        Fractional range to determine spiked pixels.
        Default (0.75, 1.25)

    Returns
    -------
    despike_image : np.ndarray
    """
    medfilt_image = ndimage.median_filter(image, size=footprint)
    spikes = image / medfilt_image
    despike_image = image.copy()
    despike_image[
        (spikes > max(srange)) | (spikes < min(srange))
    ] = medfilt_image[(spikes > max(srange)) | (spikes < min(srange))]
    return despike_image


# create a linear fit
def linear(x, a, b):
    return a + b * x


# save data to json
def save_json(data, file):
    with open(file, "w") as f:
        json.dump(data, f, indent=4)


# read data from json
def read_json(file):
    with open(file, "r") as f:
        return json.load(f)


# adjust wavelengths due to detector bias
def load_wavelengths(spectrum):
    calibration = read_json("calibration.json")
    lamda_0 = calibration["lamda_0"]
    lamda_1 = calibration["lamda_1"]

    wavelengths = []
    for i in range(len(spectrum)):
        wavelengths.append(linear(i, lamda_0, lamda_1))

    return(wavelengths)


# create a 2x bin
def bin_spectra(list):
    binned_list = []
    for i in range(int(len(list) / 2)):
        binned_list.append(list[2 * i] + list[(2 * i) + 1])
    return binned_list


#create a 2x bin
def bin_axis(spectrum, axis):
    axis = axis[0::2]
    if (len(axis) > len(spectrum)):
        axis = axis[0:len(axis) - 1]
    return axis


"""This function was provided by Dr. Sellers"""
# correct for bias, dark, flat and rotation
def load_corrections():
    with fits.open("Processing_Data/Average_Bias_20260318.fits") as hdul:
        bias = hdul[0].data
    with fits.open("Processing_Data/Average_Dark_Rate_20260318.fits") as hdul:
        darkrate = hdul[0].data
    with fits.open("Processing_Data/Processed_Gaintable_20260318.fits") as hdul:
        gain = hdul[0].data
    with fits.open("Processing_Data/Slit_Skews_20260318.fits") as hdul:
        skews = hdul[0].data
    with open("Processing_Data/beam_rotation.txt", "r") as file:
        rotation = float(file.readlines()[0])
    with open("Processing_Data/beam_edges.txt", "r") as file:
        beam_start, beam_end = file.readlines()[0].replace("\n", "").split("\t")
        beam_start = int(beam_start)
        beam_end = int(beam_end)
    return bias, darkrate, gain, skews, rotation, beam_start, beam_end


# load image and correct it
def load_data(science_file, bias, darkrate, gain):
    times = []
    with fits.open(science_file) as hdul:
        corr_data = (hdul[0].data - bias - darkrate/30) / gain
        corr_data = image_despiker(corr_data)

        header = hdul[0].header
        obs_time = get_observation_time(header)
        for t in obs_time:
            times.append(t)
    return corr_data, header, obs_time, times


"""This function was provided by Dr. Sellers"""
# rotate data
def deskew_data(corr_data, rotation, beam_start, beam_end, skews):
    derotated_data_clipped = ndimage.rotate(corr_data, rotation)[beam_start:beam_end]
    derotated_data_clipped = np.nan_to_num(derotated_data_clipped, posinf=0, neginf=0)
    deskewed_data = np.zeros(derotated_data_clipped.shape)
    for i in range(deskewed_data.shape[0]):
        deskewed_data[i] = ndimage.shift(derotated_data_clipped[i], skews[i])
    return deskewed_data


"""This function was adjusted from a function provided by Dr. Sellers"""
# rotate data and stack as a master file
def deskew_data_master(file, bias, darkrate, gain, rotation, beam_start, beam_end, skews, deskewed_data):
    with fits.open(file) as hdul:
        corr_data = (hdul[0].data - bias - (darkrate / 30)) / gain
        corr_data = image_despiker(corr_data)
    derotated_data_clipped = ndimage.rotate(corr_data, rotation)[beam_start:beam_end]
    derotated_data_clipped = np.nan_to_num(derotated_data_clipped, posinf=0, neginf=0)
    if (len(deskewed_data) == 0):
        deskewed_data = np.zeros(derotated_data_clipped.shape)
    for i in range(deskewed_data.shape[0]):
        deskewed_data[i] += ndimage.shift(derotated_data_clipped[i], skews[i])
    return deskewed_data


# load observation times
def get_observation_time(header):
    times = []
    for key in ["DATE-OBS", "DATEOBS", "TIME-OBS", "MJD-OBS", "JD", "MJD"]:
        if key in header:
            try:
                if key == "MJD-OBS" or key == "MJD":
                    obs_time = Time(header[key], format="mjd")
                elif key == "JD":
                    obs_time = Time(header[key], format="jd")
                else:
                    obs_time = Time(header[key], format="isot", scale="utc")
                times.append(obs_time.datetime)
                break
            except Exception:
                pass
    return times