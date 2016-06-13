"""
John F. Wu
2016-05-19

Here's some python code to get you started. You will find the spectra
taken by Angel Lopez-Sanchez and me at AAT/2dF+AAOmega in the directory
named aat_2015.

I'd recommend that you install a few programs and python libraries:

  * gfortran

  * linear algebra packages (blas, lapack, atlas)
    (apt-get install libblas-dev liblapack-dev libatlas-base-dev)

  * python >= 2.7
      pip        >= 1.5
      scipy      >= 0.13
      numpy      >= 1.10 (included with scipy)
      matplotlib >= 1.5  (included with scipy)
      astropy    >= 1.0 

      notebook   >= 4.2  (optional)
      ipython    >= 4.2  (optional)
      seaborn    >= 0.7  (optional)
      pandas     >= 0.12 (optional)

  * ds9 >= 7.0

  * some sort of text editor (gedit/sublime/emacs/vim)
"""


import numpy as np
from glob import glob
from astropy.wcs import WCS
from astropy.io import fits

import matplotlib.pyplot as plt
import seaborn as sns

PLOT_OPTION = 2


def read_all_data():
    '''
    Reads in the AAT spectral fits file and returns an array of the 
    spectra for each galaxy.

    Returns
    -------
    wavelengths_dict: (dict)
        a dictionary with 7 keys: 'G13', 'G14', etc.; one for each AAT
        plate. Each key corresponds to a one dimensional array of
        wavelengths of size (5037,)

    spectra_dict : (dict) 
        a dictionary with 7 keys: 'G13', 'G14', etc. Each corresponds
        to a two dimensional array of the fluxes. Each row symbolizes a
        galaxy (or fiber), and each column represents a wavelength; 
        array size is (400, 5037)
    '''

    # create list of all AAT spectra
    all_files = glob('aat_2015/G*.fits')
    all_files.sort()

    plate_ids = ['G{}'.format(n) for n in np.arange(13, 20, 1)]

    # go through each file and add their wavelengths and spectra to 
    # their respective dictionaries
    wavelengths_dict = {}
    spectra_dict = {}
    for plate_id, spectral_file in zip(plate_ids, all_files):

        spectra, header = fits.getdata(spectral_file , header=True)

        # create array of wavelengths in angstroms and add to dict
        wavelengths = np.arange(header['NAXIS1']) - header['CRPIX1']
        wavelengths = wavelengths * header['CDELT1'] + header['CRVAL1']

        wavelengths_dict[plate_id] = wavelengths

        # add spectra to dict
        spectra_dict[plate_id] = spectra

    return wavelengths_dict, spectra_dict


if __name__ == '__main__':

    wavelengths_dict, spectra_dict = read_all_data()

    # randomly choose plate G17 for now
    plate_id = 'G17'

    wavelengths = wavelengths_dict[plate_id]
    spectra     = spectra_dict[plate_id]

    # random spectra
    np.random.shuffle(spectra)

    # how many spectra to plot at once
    N = len(spectra)
    N_shown = 25

    # make a rainbow color palette
    color_palette = sns.husl_palette(N_shown, l=0.4)

    # plot every 40th galaxy
    sns.set(style='darkgrid')
    fig = plt.figure(figsize=(12, 6), dpi=150)
    ax = fig.add_subplot(111)

    if PLOT_OPTION == 1:
        
        for spectrum, color in zip(spectra[::N / N_shown], color_palette):
            ax.plot(wavelengths, spectrum, lw=0.3, color=color, alpha=0.5)

        ax.set_xlim(3750, 8950)
        ax.set_ylim(-50, 200)

    # Notice how the galaxies sort of overlap on top of each other. 
    # It might be nice to find their median flux between 6000-8000 
    # angstroms and then offset each spectrum accordingly.

    elif PLOT_OPTION == 2:
        # counter for offset; increment after each spectrum
        offset = 0  

        for spectrum, color in zip(spectra[::N / N_shown], color_palette):
            # at least one of the spectra seems to be all nans (probably 
            # guide star), so I'll skip if that's the case
            if np.isnan(spectrum).all():
                continue

            median_flux = np.nanmedian(spectrum[(wavelengths > 6000.) &
                                                (wavelengths < 8000.)])



            ax.plot(wavelengths, spectrum - median_flux + 100 * offset, 
                    lw=0.3, color = color, alpha=0.5)

            offset += 1

        ax.set_xlim(3750, 8950)
        ax.set_ylim(-50, (offset + 0.5) * 100)

    # display title and axes
    ax.set_xlabel(r'$\lambda_{\rm obs}\ [\rm \AA]$')
    ax.set_ylabel(r'$F \rm \ [arbitrary\ units]}$')
    plt.setp(ax.get_yticklabels(),visible=False)
    ax.set_title('Plate {} ({} random spectra)'.format(plate_id, offset))
    plt.tight_layout()
    plt.show()
