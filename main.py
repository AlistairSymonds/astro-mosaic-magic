import argparse
import multiprocessing
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import requests
import astropy.io.fits.header
import astropy.io.fits
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from reproject import reproject_exact, mosaicking, reproject_interp, reproject_adaptive
#from reproject.utils import reproject_blocked
from astropy.wcs import WCS
from astropy.wcs import utils as wutils
from astropy import units as u
from astropy import coordinates
from timeit import default_timer as timer
import pprint
import json


def calc_sky_area(header):
    return wutils.proj_plane_pixel_area(WCS(header)) * header['IMAGEH'] * header['IMAGEW']

def calc_px_scale(pixel_size_um, focal_length_mm):
    return (pixel_size_um / focal_length_mm) * 206.265

class MosaicMagic:


    def __init__(self, astrometry_url):
        self.as_url = astrometry_url
        login_resp = requests.post(astrometry_url+"/api/login")
        login_data = json.loads(login_resp.text)

        #set up astrometry instance
        self.session_id = login_data["session"]
        print("Got astrometry session id: " + self.session_id)



    def solve_image(self, img, force_solve=False):
        print(img)

    def process_single_image(self, fits_paths: Path):
        print("Starting to process image at "  + str(fits_paths.absolute()))
        hdul = astropy.io.fits.open(fits_paths.absolute())
        hdul.info()
        pprint.pprint(hdul[0].header)

        px_scale = calc_px_scale(hdul[0].header['XPIXSZ'], hdul[0].header['FOCALLEN'] )
        
        blind_solve = False
        ra = 0
        dec = 0

        if hdul[0].header['RA'] is not None:
            ra = hdul[0].header['RA']
        else:
            blind_solve = True

        if hdul[0].header['DEC'] is not None:
            dec = hdul[0].header['DEC']
        else:
            blind_solve = True
        
        print(f"Estimated pixel scale from fits headers is {px_scale}")
        if not blind_solve:
            print(f"Extracted following RA: {ra} and DEC: {dec}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--imgs_dir", help="Path to folder containing fits for registration", default="imgs")
    ap.add_argument("--arc_sec_per_px", help="An override for final arc sec/px, if no value is given the smallest "
                                             "values in all images will be used")
    ap.add_argument("--astrometry_port", default="8187")
    ap.add_argument("--swap_dir", help="Directory memory mapped ararys will be saved during stacking. "
                                       "A fast SSD with lots of space will increase speed, however having enough space is more important",
                    default="mosaic_swap")

    args = ap.parse_args()

    px_scale = None
    if args.arc_sec_per_px is not None:
        px_scale = float(args.arc_sec_per_px) * u.arcsec


    # check args
    img_dir_path = None
    if args.imgs_dir is not None:
        img_dir_path = Path(args.imgs_dir)
        if img_dir_path.exists() == False or img_dir_path.exists() == False:
            assert NotADirectoryError("Couldn't find directory: "+ str(args.imgs_dir))
   
    fits_for_processing = []
    for fits_path in img_dir_path.glob("*.fits"):
        fits_for_processing.append(fits_path)

    if len(fits_for_processing) == 0:
        print("Couldn't find any .fits files under path " + str(img_dir_path.absolute()))
        exit()

    pprint.pprint(fits_for_processing)
    print("Found " + str(len(fits_for_processing)) + " .fits files to process and mosiac")

    base_url = "http://127.0.0.1:"+str(args.astrometry_port)

    mm = MosaicMagic(base_url)
    

    solve_pool = ProcessPoolExecutor()
    #processing_futures = []
    #for fits in fits_for_processing:
    #    processing_futures.append(solve_pool.submit(process_single_image, fits))
    
    mm.process_single_image(fits_for_processing[0])

    print("Starting to wait")
    #concurrent.futures.wait(processing_futures)
    print("All platesolving done")



if __name__ == "__main__":
    multiprocessing.freeze_support()
    Image.MAX_IMAGE_PIXELS = None #keep PIL happy when opening stupidly 3x drizzled files :P
    main()