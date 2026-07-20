import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import tifffile

from astropy.io import fits
from astropy.time import Time
from photutils.centroids import centroid_2dg
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry


# measure flux of single exposure, adjust aperature and annuli as needed
def measure_bright_star_flux(
    fits_file,
    aperture_radius=25.0,
    annulus_r_in=50.0,
    annulus_r_out=100.0,
    cutout_half_size=25,
    position=[]
):
    # load file
    with fits.open(fits_file) as hdul:
        # load corrections
        bias, dark, flat = load_corrections()
        data = (hdul[0].data - bias - dark)/flat
        header = hdul[0].header

    # initialize positions of stars
    x_pos = [position[0][0], position[0][1]]
    y_pos = [position[1][0], position[1][1]]
    temp = data[y_pos[0]:y_pos[1], x_pos[0]:x_pos[1]]

    # find brightest pixel
    y_peak, x_peak = np.unravel_index(np.argmax(temp), temp.shape)
    y_peak = y_peak + y_pos[0]
    x_peak = x_peak + x_pos[0]

    # ake a cutout around the brightest pixel
    y_min = max(0, y_peak - cutout_half_size)
    y_max = min(data.shape[0], y_peak + cutout_half_size + 1)
    x_min = max(0, x_peak - cutout_half_size)
    x_max = min(data.shape[1], x_peak + cutout_half_size + 1)

    cutout = data[y_min:y_max, x_min:x_max]

    # estimate center in the cutout
    try:
        x_c_local, y_c_local = centroid_2dg(cutout)
        x_center = x_min + x_c_local
        y_center = y_min + y_c_local
    except Exception:
        # use brightest pixel if centroid fit fails
        x_center = float(x_peak)
        y_center = float(y_peak)

    position = [(x_center, y_center)]

    # set aperture and annulus
    aperture = CircularAperture(position, r=aperture_radius)
    annulus = CircularAnnulus(position, r_in=annulus_r_in, r_out=annulus_r_out)

    # sum aperture
    ap_table = aperture_photometry(data, aperture)
    raw_aperture_sum = float(ap_table["aperture_sum"][0])

    # calculate background from annulus
    annulus_masks = annulus.to_mask(method="center")
    annulus_data = annulus_masks[0].multiply(data)
    annulus_pixels = annulus_data[annulus_masks[0].data > 0]

    # exclude hot pixels
    annulus_pixels = annulus_pixels[np.isfinite(annulus_pixels)]

    # calculate median background
    background_median = float(np.median(annulus_pixels))

    # calculate flux from background
    aperture_area = float(aperture.area)
    background_in_aperture = background_median * aperture_area

    # remove background
    net_flux = raw_aperture_sum - background_in_aperture

    # read observation time
    obs_time = None
    for key in ["DATE-OBS", "DATEOBS", "TIME-OBS", "MJD-OBS", "JD", "MJD"]:
        if key in header:
            try:
                if key == "MJD-OBS" or key == "MJD":
                    obs_time = Time(header[key], format="mjd")
                elif key == "JD":
                    obs_time = Time(header[key], format="jd")
                else:
                    obs_time = Time(header[key], format="isot", scale="utc")
                break
            except Exception:
                pass

    # create object to plot
    return {
        "file": os.path.basename(fits_file),
        "x_center": x_center,
        "y_center": y_center,
        "raw_aperture_sum": raw_aperture_sum,
        "background_median_per_pixel": background_median,
        "background_in_aperture": background_in_aperture,
        "net_flux": net_flux,
        "obs_time": obs_time,
    }


# create a list for better handling
def create_list(star_data):
    data = []
    for measurement in star_data:
        data.append(measurement["net_flux"])
    return data

# load bias, dark and flat corrections
def load_corrections():
    with fits.open("master_bias.fits") as hdul:
        bias = hdul[0].data
    with fits.open("master_dark_bias_corrected.fits") as hdul:
        dark = hdul[0].data

    flat = tifffile.imread("master_flat.tif").astype(float)
    flat /= flat.mean()

    return bias, dark, flat

# correct dimming by normalizing at constant star
def correct_dimming(variable, fixed):
    variable = np.asarray(create_list(variable))
    fixed = np.asarray(create_list(fixed))
    return variable / fixed

# process whole object folder
def process_folder(
    folder="AD_Leo",
    aperture_radius=8.0,
    annulus_r_in=12.0,
    annulus_r_out=18.0,
    cutout_half_size=15,
    positions=[]
):
    # load files
    fits_files = sorted(glob.glob(os.path.join(folder, "*.fits")))

    results = []

    counter = 0
    for pos in positions:
        counter += 1
        pos_results = []
        # select necessary data
        for fits_file in fits_files[0:1600]:
            try:
                # measure flux of exposure
                result = measure_bright_star_flux(
                    fits_file,
                    aperture_radius=aperture_radius,
                    annulus_r_in=annulus_r_in,
                    annulus_r_out=annulus_r_out,
                    cutout_half_size=cutout_half_size,
                    position=pos
                )
                pos_results.append(result)
                print(f"Run {counter} : : {result['file']}: net_flux = {result['net_flux']:.3f}")
            except Exception as e:
                print(f"Skipping {os.path.basename(fits_file)} because of error: {e}")
        results.append(pos_results)

    # if all observations have a time, use time on x axis
    all_have_time = all(r["obs_time"] is not None for r in results[0])

    # run for AD_Leo and reference star
    ad_leo = results[0]
    s1 = results[1]
    stars = [ad_leo, s1]

    x_array = []
    flux_array = []

    # insert x axis data
    if all_have_time:
        for star in stars:
            star.sort(key=lambda r: r["obs_time"].jd)
            x_array.append([r["obs_time"].datetime for r in star])
            x_label = "Time [UTC]"
            use_dates = True
    else:
        for star in stars:
            x_array.append(np.arange(len(star)))
            x_label = "File index"
            use_dates = False

    # add fluxes for y axis
    for star in stars:
        fluxes = [r["net_flux"] for r in star]
        flux_array.append(fluxes)

    # plot brightness curves
    plt.figure(figsize=(8, 5))
    for i in range(2):
        plt.scatter(x_array[i], flux_array[i], marker=".", linestyle="-")
    plt.xlabel(x_label)
    plt.ylabel("Intensity")
    plt.title("Change in intensity")
    plt.legend(["AD Leonis", "TYC 1423-165-1"])

    # set x axis to dates if possible
    if use_dates:
        plt.gcf().autofmt_xdate()
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    plt.tight_layout()
    plt.show()

    # plot normalized intensity
    flux = correct_dimming(stars[0], stars[1])

    plt.figure(figsize=(8, 5))
    plt.scatter(x_array[0], flux, marker=".", linestyle="-")
    plt.xlabel(x_label)
    plt.ylabel("Intensity")
    plt.title("Normalized Change in AD Leonis Intensity")

    # set x axis to dates if possible
    if use_dates:
        plt.gcf().autofmt_xdate()
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    plt.tight_layout()
    plt.show()
    return results


# calculate intensity for whole data set
if __name__ == "__main__":
    # position of AD_Leo and reference star
    positions = [[[1500, 1900], [800, 1200]], [[1550, 2050], [1050, 1500]]]
    results = process_folder(
        folder="AD_Leo",
        aperture_radius=25.0,
        annulus_r_in=50.0,
        annulus_r_out=100.0,
        cutout_half_size=25,
        positions=positions
    )