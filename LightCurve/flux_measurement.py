import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import tifffile
from astropy.io import fits
from photutils.centroids import centroid_2dg
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry


# measure flux of single exposure, adjust aperature and annuli as needed
def measure_bright_star_flux(
    fits_file,
    aperture_radius=25.0,
    annulus_r_in=50.0,
    annulus_r_out=100.0,
    cutout_half_size=25,
    positions=[],
    show_plot=True
):
    # load file
    with fits.open(fits_file) as hdul:
        # load and apply corrections
        bias, dark, flat = load_corrections()
        data = (hdul[0].data - bias - dark) / flat

    flux_list = []
    results = []

    # initialize positions of stars
    for pos in positions:
        x_pos = [pos[0][0], pos[0][1]]
        y_pos = [pos[1][0], pos[1][1]]
        temp = data[y_pos[0]:y_pos[1], x_pos[0]:x_pos[1]]

        # find brightest pixel
        y_peak, x_peak = np.unravel_index(np.argmax(temp), temp.shape)
        y_peak = y_peak + y_pos[0]
        x_peak = x_peak + x_pos[0]

        # make a cutout around the brightest pixel
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
        raw_aperture_sum = ap_table["aperture_sum"][0]

        # calculate background from annulus
        annulus_masks = annulus.to_mask(method="center")
        annulus_data = annulus_masks[0].multiply(data)
        annulus_pixels = annulus_data[annulus_masks[0].data > 0]

        # exclude hot pixels
        annulus_pixels = annulus_pixels[np.isfinite(annulus_pixels)]

        # calculate median background
        background_median = np.median(annulus_pixels)

        # calculate flux from background
        aperture_area = aperture.area
        background_in_aperture = background_median * aperture_area

        # remove background
        net_flux = raw_aperture_sum - background_in_aperture

        flux_list.append(net_flux)

        # create object to plot
        result = {
            "x_center": x_center,
            "y_center": y_center,
            "raw_aperture_sum": raw_aperture_sum,
            "background_median_per_pixel": background_median,
            "background_in_aperture": background_in_aperture,
            "net_flux": net_flux,
            "x_pos": x_pos,
            "y_pos": y_pos,
            "aperture": aperture,
            "annulus": annulus
        }
        results.append(result)

    # plot a single exposure as demonstration of process
    if show_plot:
        plt.figure(figsize=(8, 5))
        vmin, vmax = np.percentile(data, [5, 99.5])
        plt.imshow(data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)

        plt.title("AD Leonis flux measurement")
        ax = plt.gca()
        for res in results:
            # plot annuli and rectangle
            width = res["x_pos"][1] - res["x_pos"][0]
            height = res["y_pos"][1] - res["y_pos"][0]
            rect = Rectangle((res["x_pos"][0], res["y_pos"][0]), width, height, edgecolor="r", facecolor="none", linewidth=1)
            ax.add_patch(rect)
            res["aperture"].plot(color="red", lw=1)
            res["annulus"].plot(color="blue", lw=1)
        plt.xlabel("x [px]")
        plt.ylabel("y [px]")
        plt.tight_layout()
        plt.show()

    return results


# load bias, master and flat correction
def load_corrections():
    with fits.open("master_bias.fits") as hdul:
        bias = hdul[0].data
    with fits.open("master_dark_bias_corrected.fits") as hdul:
        dark = hdul[0].data

    flat = tifffile.imread("master_flat.tif").astype(float)
    flat /= flat.mean()

    return bias, dark, flat

# create plot of single exposure
if __name__ == "__main__":
    fits_file = "AD_Leo/AD Leonis_2026-04-02_00-29-44__-25.00_1.00s_0629.fits"

    # adjust positions if exposure needs it, these work for AD_Leo
    positions = [[[1500, 1900], [800, 1200]], [[1550, 2050], [1050, 1500]]]

    # run measurement
    result = measure_bright_star_flux(
        fits_file,
        aperture_radius=25.0,
        annulus_r_in=50.0,
        annulus_r_out=100.0,
        cutout_half_size=25,
        positions=positions,
        show_plot=True
    )
