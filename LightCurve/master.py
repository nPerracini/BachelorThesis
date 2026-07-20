from astropy.io import fits
import numpy as np
import glob
import os


# create master frames for brightness curve
def create_master_frame(folder):

    #load files
    files = glob.glob(os.path.join(folder, "*.fits"))

    data_stack = []

    # read files
    for file in files:
        with fits.open(file) as hdul:
            data_stack.append(hdul[0].data.astype(float))

    data_stack = np.array(data_stack)

    #create master frame
    master_frame = np.median(data_stack, axis=0)

    return master_frame

# create master bias
def create_bias():
    master_bias = create_master_frame("Bias")

    # save master bias
    fits.writeto("master_bias.fits", master_bias, overwrite=True)

# create bias corrected master dark
def create_dark():
    # load master bias
    master_bias = fits.getdata("master_bias.fits")

    dark_files = glob.glob("Dark/*.fits")

    dark_stack = []

    for file in dark_files:
        dark = fits.getdata(file).astype(float)

        # remove bias from each dark
        dark -= master_bias

        dark_stack.append(dark)

    master_dark = np.median(dark_stack, axis=0)

    # save bias corrected master dark
    fits.writeto("master_dark_bias_corrected.fits",
                 master_dark,
                 overwrite=True)

# run needed functions
create_bias()
create_dark()