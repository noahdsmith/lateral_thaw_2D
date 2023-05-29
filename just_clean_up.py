import xarray as xa
import argparse
import numpy as np
import time, os, glob

from global_lat_thaw_parm_mod import output_path, prefix_prefix

parser = argparse.ArgumentParser()
parser.add_argument('step', type=int, help='Longitude step')
args = parser.parse_args()
step = args.step

lons = np.arange(-180,180,args.step)

print('Warning - will delete stuff!')

print(f'Looking in folder : {output_path}')

if input("Proceed? (T/F) : ") == 'T':

    print('Cleaning up')
    for lon in lons:
        try:
            os.remove(f'{output_path}{prefix_prefix}_{lon:+04}_{lon+step:+04}year.nc')
        except OSError:
            print(f"Couldn't clean up {output_path}{prefix_prefix}_{lon:+04}_{lon+step:+04}year.nc")
        file_list = glob.glob(f'{output_path}{prefix_prefix}_{lon:+04}_{lon+step:+04}month*.nc')
        for filename in file_list:
            try:
                os.remove(filename)
            except OSError:
                print(f"Couldn't clean up {filename}")

    print('Done :)')