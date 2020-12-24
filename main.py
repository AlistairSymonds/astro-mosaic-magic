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
from astroquery.astrometry_net import AstrometryNet


def calc_sky_area(header):
    return wutils.proj_plane_pixel_area(WCS(header)) * header['IMAGEH'] * header['IMAGEW']

def calc_px_scale(pixel_size_um, focal_length_mm):
    return (pixel_size_um / focal_length_mm) * 206.265

class MosaicMagic:


    def __init__(self, api_key, astrometry_url):
        self.as_url = astrometry_url
        self.ast = AstrometryNet()
        self.ast.api_key = api_key

        

        #login_resp = requests.post(self.as_url+"/api/login")
        #login_data = json.loads(login_resp.text)

        #set up astrometry instance
        #self.session_id = login_data["session"]
        #print("Got astrometry session id: " + self.session_id)

    def upload_img(self, fits_paths: Path, img, arcsecperpx, scale_error=100):

        headers = {'Content-type': 'multipart/form-data'} 
        
        solve_meta = {}
        
        #header["Content-Type"] = "multipart/form-data"
        solve_meta["session"] = self.session_id
        solve_meta["allow_commercial_use"] = "d"
        solve_meta["allow_modifications"] = "d"
        solve_meta["publicly_visible"] = "y"
        solve_meta["scale_units"] = "arcsecperpix"
        solve_meta["scale_type"] = "ev"
        solve_meta["scale_est"] = arcsecperpx
        solve_meta["scale_error"] = scale_error

        header_json = json.dumps(solve_meta)
        print(header_json)

        
        files = {'file': open(fits_paths, 'rb')}
        print("uploading")
        test = requests.get(self.as_url+"/api/upload/", data=solve_meta)

        print(test.text)


    def solve_image(self, fits_paths: Path, img, force_solve=False):
        print(img)

        px_scale = calc_px_scale(img.header['XPIXSZ'], img.header['FOCALLEN'] )
        
        blind_solve = False
        ra = 0
        dec = 0

        if img.header['RA'] is not None:
            ra = img.header['RA']
        else:
            blind_solve = True

        if img.header['DEC'] is not None:
            dec = img.header['DEC']
        else:
            blind_solve = True
        
        print(f"Estimated pixel scale from fits headers is {px_scale}")
        if not blind_solve:
            print(f"Extracted following RA: {ra} and DEC: {dec}")

        existing_wcs = WCS(img.header)
        if not existing_wcs.has_celestial and not force_solve:
            #self.upload_img(fits_paths, img, px_scale)
            try_again = True
            submission_id = None

            while try_again:
                try:
                    if not submission_id:
                        wcs_header = self.ast.solve_from_image(str(fits_paths),
                                                        submission_id=submission_id,
                                                        solve_timeout=1200)
                    else:
                        wcs_header = self.ast.monitor_submission(submission_id,
                                                            solve_timeout=600)
                except TimeoutError as e:
                    submission_id = e.args[1]
                else:
                    # got a result, so terminate
                    try_again = False

            
            w = WCS(wcs_header)
            solved_w_header = w.to_fits(relax=True)
            print(solved_w_header[0])
            for header_key in solved_w_header[0].header:
                print(f"Adding [{header_key}] = {solved_w_header[0].header[header_key]}")
                if header_key in img.header:
                    print("Print key already exists in dest file, skipping")
                else:
                    img.header[header_key] = solved_w_header[0].header[header_key]
        else:
            print("Already detected valid WCS header data in " + str(fits_paths.absolute()) + " - skipping solve!")

        pprint.pprint(img.header)

        return img

    def process_single_image(self, fits_paths: Path):
        print("Starting to process image at "  + str(fits_paths.absolute()))
        hdul = astropy.io.fits.open(fits_paths.absolute(), mode='update')
        pprint.pprint(hdul[0].header)

        hdul[0] = self.solve_image(fits_paths, hdul[0])

        hdul.flush()
        hdul.close()
    

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--imgs_dir", help="Path to folder containing fits for registration", default="imgs")
    ap.add_argument("--arc_sec_per_px", help="An override for final arc sec/px, if no value is given the smallest "
                                             "values in all images will be used")
    ap.add_argument("--astrometry_port", default="8187")
    ap.add_argument("--swap_dir", help="Directory memory mapped ararys will be saved during stacking. "
                                       "A fast SSD with lots of space will increase speed, however having enough space is more important",
                    default="mosaic_swap")
    ap.add_argument("--astrometry_api_key")

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

    mm = MosaicMagic(args.astrometry_api_key, base_url)
    

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