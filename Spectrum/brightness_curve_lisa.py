import os

from functions import load_corrections
from functions import load_data
from functions import deskew_data

import matplotlib.pyplot as plt
import matplotlib.dates as mdates


final_data = []
times = []

# create brightness curve from provided folder
if __name__ == "__main__":
    folder = "Data/AD_Leonis_new"

    # sort files by date
    files = sorted(
        os.listdir(folder),
        key=lambda f: os.path.getctime(os.path.join(folder, f))
    )

    # process each file
    for spectra in files:
        science_file = os.path.join(folder, spectra)

        bias, darkrate, gain, skews, rotation, beam_start, beam_end = load_corrections()

        corr_data, header, obs_time, t = load_data(science_file, bias, darkrate, gain)

        deskewed_data = deskew_data(corr_data, rotation, beam_start, beam_end, skews)

        times.append(t)
        final_data.append(deskewed_data.sum(axis=0)[15:-15].sum() / len(deskewed_data[0]))
        print("Finished for " + spectra)

    # plot brightness curve
    plt.figure(figsize=(8, 5))
    plt.plot(times, final_data, linewidth=1)
    plt.gcf().autofmt_xdate()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.xlabel("Time [UTC]")
    plt.ylabel("Intensity")
    plt.title("Change in AD Leonis intensity")
    plt.show()
