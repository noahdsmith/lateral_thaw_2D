import xarray as xa
import argparse
import numpy as np
import time, os, glob

from global_lat_thaw_parm_mod import output_path, prefix_prefix, output_profiles

parser = argparse.ArgumentParser()
parser.add_argument('n_tasks', type=int, help='Number of tasks')
args = parser.parse_args()
n_tasks = args.n_tasks

test_verbose = False # just print what would happen

def concat_output(output_suffix, concat_type):
    print(f"Opening & concatenating datasets for {output_suffix}")
    start = time.time()
    try:
        if test_verbose:
            print(f'Open and concatenate {output_path}{prefix_prefix}_???{output_suffix}.nc')
            print(f'Output as {output_path}{prefix_prefix}{output_suffix}.concat.nc')
        else:
            output_data = xa.open_mfdataset(f'{output_path}{prefix_prefix}_???{output_suffix}.nc',concat_dim = 'land')

            end = time.time()
            print(f'Opened & concatenated in {end-start}s')

            print('Saving output')
            start = time.time()

            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in output_data.data_vars}
            output_data.to_netcdf(output_path + prefix_prefix+f"{output_suffix}.concat.nc", encoding=encoding)

            end = time.time()
            print(f'Saved in {end-start}s')

        print('Cleaning up')
        file_list = glob.glob(f'{output_path}{prefix_prefix}_??{output_suffix}.nc')
        if test_verbose:
            print('Remove ', file_list)
        else:
            for filename in file_list:
                try:
                    os.remove(filename)
                except OSError:
                    print(f"Couldn't clean up {filename}")
        if concat_type == 'daily':
            file_list = glob.glob(f"{output_path}{prefix_prefix}_???{output_suffix}month*.nc")
            if test_verbose:
                print('Remove ', file_list)
            else:
                for filename in file_list:
                    try:
                        os.remove(filename)
                    except OSError:
                        print(f"Couldn't clean up {filename}")

    except:
        print(f"A problem occured for {output_suffix}")

for output in output_profiles:
    concat_output(f"{output['freq']}_{output['vars']}",output['freq'])

concat_output('.setup_dump','setup_dump')

print('Done :)')