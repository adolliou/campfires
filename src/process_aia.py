#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import astropy
import astropy.units as u
from astropy.io import fits
from astropy.io.fits import CompImageHDU
import sunpy.map
from sunpy.net import Fido, attrs as a
from aiapy.calibrate import register, update_pointing
from aiapy.calibrate.util import get_pointing_table
from pathlib import Path
import os

#aia_data_path  = r'/home/francisv/CLOSE-UP/new_base_code_20250911/test_run_20200530/data/aia/171/'
#start_time     = r'2020-05-30T14:54:00'
#end_time       = r'2020-05-30T14:58:20'
#aia_wavelength = 171  # can be 171, 94, 131, 193, 211, 335


def process_aia(aia_data_path, start_time, end_time, wavelength):
    # use sunpy tool Fido to find and download level 1 AIA data
    q = Fido.search(a.Time(start_time, end_time),
                    a.Instrument('AIA'),
                    a.Wavelength(wavemin = wavelength*u.angstrom, wavemax = wavelength*u.angstrom))
    
    output_dir = Path(aia_data_path)
    output_dir.mkdir(exist_ok=True)
   #downloaded_files = Fido.fetch(q, path = aia_data_path)
    downloaded_files = Fido.fetch(q, path = aia_data_path, site = "NSO")
    print(f"Downloaded {len(downloaded_files)} FITS files")

   #correction_table = get_correction_table("JSOC")  # needed in case you want to correct for instrument degradation

    for i, f in enumerate(downloaded_files, start=1):
        try:
            aia_map = sunpy.map.Map(f)

            with fits.open(f) as hdul:
                if (any(isinstance(hdu, CompImageHDU) for hdu in hdul)):
                    header = fits.getheader(f, 1)  # for compressed fits files
                else:
                    header = fits.getheader(f)  # for uncompressed fits files

            print("Obtained header")

            if ("OSCNMEAN" in header):
                del header["OSCNMEAN"]
            if ("OSCNRMS" in header):
                del header["OSCNRMS"]

            # prep the level 1 AIA file to level 1.5
            pointing_table = get_pointing_table("JSOC", time_range=(aia_map.date - 12 * u.h, aia_map.date + 12 * u.h))
            aia_map_updated_pointing = update_pointing(aia_map, pointing_table=pointing_table)
            aia_prep_map = register(aia_map_updated_pointing)  # The image in aia_prep_map is now a level 1.5 data product
           #aia_corrected_map = correct_degradation(aia_prep_map, correction_table = correction_table)  # needed in case you want to correct for instrument degradation
            
            base = os.path.basename(f)
            out = os.path.join(aia_data_path, base.replace("lev1.fits", "lev1p5.fits"))
            aia_prep_map.save(out, overwrite=True)
            print(f"[{i}/{len(downloaded_files)}] Saved: {out}")
            
        except Exception as e:
            print(f"[{i}/{len(downloaded_files)}] Failed to process {f}: {e}")
        print("All files processed and saved to:", output_dir)
