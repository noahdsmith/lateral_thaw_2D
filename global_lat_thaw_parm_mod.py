#----------------------------------
# Imports
#----------------------------------
print('Importing modules and data')

code_test = False

print('global_lat_thaw_parm_mod last updated 30/03/2023 17:22')

###
# Notes
# --------
# Now that interpolating temperatures, could totally clean up alot of the driving options?
###

import xarray as xa
import numpy as np
import datetime
import pandas as pd
#import calendar
import time
import pathlib
import pickle # needed here?
from tqdm import tqdm
from scipy import interpolate

#----------------------------------
# Parameters
#----------------------------------

# Set for each run
#------------------------
default_parm_options = dict(
        option = 'offset', # twotilesat or pardry, used to be a highanddry too.
duo = True, # check this works elsewhere!
hist_or_sspno = 'historical', # for future runs [ ssp585, ssp370, ssp126]
surf_daily = False, # check this works elsewhere! # note that the most recent duo is monthly only.
crop = 'land', # 'lon'(broken), 'latlon' or 'land'
minlon = 25.0, maxlon = 26.0, minlat = 68.0, maxlat = 70.0, # if cropping by latlon
soil_prop_to_use = 'peat', # organic or peat, used to be a dynamic dump input too
mineral_layer = True, # if using offset and excess ice
mineral_depth_bs = 0.75, # was 1.5
xice_content = 0.25, # was 0.3
xice_depth_bs = 0.75, # Usually set to top of mineral layer.
high_bit_elevation = 0.6, # Palsa height (m) Warning, things might break if changed!
slope_x = 1.2, # The horizontal length of the slope (m)
palsa_width = 6.0, # Width of flat top of palsa (m) , used if grid_option == 'double_sided'.
temp_adjust = False, # This should only be used for sat_unsat 2010 and uses 2090 to 2100 MAT for ssp5.8
excess_ice = True,
grid_option = 'dense', # was false, might be a good idea to use with excess ice?
depth_beyond_7 = 0, # extra 1 m thick layers added below the usual 7 m depth
T_shift = 0.0, # Doesn't work for all driving. Just to test.
param_option = 'default',
y_threshold_lower = 0.15, # Surface elevation threshold for x_of_y_threshold diagnostic.
y_threshold_upper = 0.50,
domain_x_threshold = -0.4,
shifting_domain = True # Shifts model domain based on x_of_y_threshold to avoid edge effects. Vulnerable to the code changing.
    ) # sim_parms

class thaw_sim(object): # At this stage, mainly because it makes the syntax easier
    def __init__(self, initial_data):
        for key in initial_data:
            setattr(self, key, initial_data[key])

if not code_test:
    # Terrible, but do a classier approach some other time.
    parm_pickle_file = open('sel_parm_options.pk', 'rb')
    parm_options_sel = pickle.load(parm_pickle_file)
    for key in parm_options_sel:
        print(f'Overwriting {key} with {parm_options_sel[key]}')
        default_parm_options[key] = parm_options_sel[key]
sp = thaw_sim(default_parm_options)
if not code_test:
    parm_pickle_file.close()

if sp.shifting_domain and (sp.grid_option not in ['dense', 'dense_longer']):
    print(f'Warning! grid option {sp.grid_option} chosen, make sure this has constant x spacing')

excess_ice_code = 'original' # 'original' is reccomended. 'multiarray' is also available, but is weirdly slower, and may have bugs.
interpolate_sthu = True
first_order_interp = True # More smearing but avoids instabilities, doesn't work for multiarray.
time_to_go = 'print_daily' # Options are 'tqdm', 'print_daily' or 'silent'.
correct_thaw = True # Correct for xice thaw numerics

#output_path = '/work/scratch-pw/nds211/batch_run/10_Iskoras_offset_1850_1910/'
if code_test:
    output_path = 'test_output/'
else:
    output_path =  f'/work/scratch-pw3/nds211/batch_run/13_xice_shift_options_2/{sp.hist_or_sspno}/' # f'/gws/nopw/j04/jules/noahsmith/batch_run/12_xice_{sp.hist_or_sspno}_global/' # sp.grid_option
print(f'Output_path = {output_path}')
prefix_prefix = f'13_global_xice_{sp.hist_or_sspno}_{sp.grid_option}_'
output_profiles =  [{'freq':'monthly','vars':'full'}]# # min_xice, full
                    #{'freq':'monthly','vars':'ends'}]
                   # {'freq':'yearly','vars':'reduced_year'}] # set here
avg_t_soil = True

if sp.surf_daily and sp.option != 'pardry':
    print('Warning: since surf_daily, usually you need pardry?') # noheatflux edges too?
    
if sp.surf_daily and sp.duo:
    print('Warning: there was a problem with the daily duo data being monthly.') # noheatflux edges too?

if sp.high_bit_elevation % 0.1 != 0.0:
    print('Warning, at the moment code needs high_bit_elevation % 0.1 == 0.0')

# Other settings
#-----------------------
no_heat_flux_edges = True # MUST BE TRUE FOR LAT THAW CALC TO WORK

# Temperature change
#------------------------
dpsidt = 114.3 # Rate of change of ice potential with temperature
               # RHO_ICE*LF/ZERODEGC*1/(RHO_WATER*G) (m/K)
zerodegc = 273.15
rho_water = 1000.0 # density of pure water (kg/m3)
hcapi = 2100.0 # Specific heat capacity of ice (J/kg/K)
hcapw= 4180.0 # Specific heat capacity of water (J/kg/K)
lf = 0.334e6 # latent heat of fusion of water at 0degc (J kg-1)

# Conductivity
#------------------------
hcair = 0.025 # Thermal conductivity of air (W/m/K)
hcice = 2.24  # Thermal conductivity of ice (W/m/K)
hcwat = 0.56  # Thermal conductivity of liquid water (W/m/K)

# Numerics
#------------------------
tiny_0 = np.finfo(float).eps #tiny(0.0)
small_value = np.finfo(float).eps # epsilon(0.0)
    # Not sure about these two
    
# Time
#------------------------
timestep = 3600 # Default = 3600 seconds (hourly). Saving will break for a timestep >= daily

# Space
#------------------------
thickness = 0.4 # probably! # dz[n] ! Thicknesses of the cell (m).

if sp.option == 'twotilesat':
    nx, ny = 2,20 + sp.depth_beyond_7
    
elif sp.option == 'offset':
    std_x_spacings = {
        'standard' : np.array([0.4]*8 + [0.2] * 8 + [0.4] * 8),
        'dense' : np.array([0.2]*16 + [0.2] * 8 + [0.2] * 16),
        'dense_longer' : np.array([0.2]*10 + [0.2]*16 + [0.2] * 8 + [0.2] * 16 + [0.2]*10),
        'half_dense_longer': np.array([0.2]*32 + [0.2] * 8 + [0.4] * 16),
        'double_sided': np.array([0.4] * 8 + [0.2]*40 + [0.4] * 8)
    }[sp.grid_option]

    std_y_spacings = np.array([0.1]*int(sp.high_bit_elevation/0.1) + [0.1]*6 + [0.2]*5 + [0.3]*5 + [0.8] * 6 + [1.0]*sp.depth_beyond_7) # have turned to array, hopefully doesn't break anything!
    
    def thickness_to_positions(thicknesses = None, center = False):
        """ For converting thicknesses of layers to midpoint positions. """
        positions = thicknesses/2
        positions[1:] += np.cumsum(thicknesses[:-1])
        if center:
            positions = positions - sum(thicknesses)/2
        return positions
    
    nx, ny = len(std_x_spacings), len(std_y_spacings)
    
    std_x_positions = thickness_to_positions(std_x_spacings,center = True)
    std_y_positions = sp.high_bit_elevation-thickness_to_positions(std_y_spacings)
    tot_depth = std_y_positions[-1] - std_y_spacings[-1]/2.0
    xx_temp_o, yy_temp_o = np.meshgrid(std_x_positions, std_y_positions, indexing='ij')
    dx_xymatrix, dy_xymatrix = np.meshgrid(std_x_spacings, std_y_spacings, indexing='ij')
 
    slope_x = sp.slope_x # horizontal slope length 
    slope_y = sp.high_bit_elevation # height of palsa
    surface_thickness = std_y_spacings[0] # Edge to edge thickness of surface layer
    
    def get_surface_height(x, slope_x, slope_y, grid_option, palsa_width):
        if grid_option != 'double_sided':
            if x < -slope_x/2.0:
                h = slope_y
            elif x > slope_x/2.0:
                h = 0.0
            else:
                h = slope_y * (0.5 - (x/slope_x))
        else:
            if x > -palsa_width/2.0 and x < palsa_width/2.0:
                h = slope_y
            elif x < -palsa_width/2.0 - slope_x or x > palsa_width/2.0 + slope_x:
                h = 0.0
            else:
                h = slope_y * (0.5 - ((np.sqrt(x**2)-(palsa_width+slope_x)/2)/slope_x))
        return h
    if sp.grid_option != 'double_sided':
        left_slope = np.arange(nx)[std_x_positions>(-slope_x/2.0)].min()
        right_slope = np.arange(nx)[std_x_positions>slope_x/2.0].min() # Careful! Old options like pardry still overwrite these.
    else: # For double sided, these are the indices on the positive side of the palsa
        left_slope = np.arange(nx)[std_x_positions>(sp.palsa_width/2.0)].min()
        right_slope = np.arange(nx)[std_x_positions>(sp.palsa_width/2.0 + slope_x)].min() # Careful! Old options like pardry still overwrite these.

    depth_bs = yy_temp_o.copy() # depth below surface
    xy_mask = np.ones_like(yy_temp_o)
    temp_y_indices = np.indices(np.shape(xx_temp_o))[1]
    for i, j in zip(np.indices(np.shape(xx_temp_o))[0].flatten(), np.indices(np.shape(xx_temp_o))[1].flatten()):
        if yy_temp_o[i,j] + dy_xymatrix[i,j]/2 > get_surface_height(xx_temp_o[i,j],slope_x, slope_y, sp.grid_option, sp.palsa_width):
            depth_bs[i,j] = np.nan # presumably could just times the mask by these
            xy_mask[i,j] = np.nan
            dx_xymatrix[i,j], dy_xymatrix[i,j] = np.nan, np.nan
            temp_y_indices[i,j] = 3e6
    surf_y_index = np.nanmin(temp_y_indices, axis = 1)
    std_layer_tops = std_y_positions + np.array(std_y_spacings)/2
    std_layer_bottoms = std_y_positions - np.array(std_y_spacings)/2
    #surface_heights = np.array([get_surface_height(x, slope_x, slope_y, sp.grid_option, sp.palsa_width) for x in std_x_positions])
    surface_heights = np.array([std_layer_tops[surf_y_index[i]] for i in range(nx)]) # for some reason different to the above! maybe an inequality.
    depth_bs = np.transpose(np.transpose(depth_bs) - surface_heights) # transposing to broadcast along x
    
    if sp.excess_ice:
        # This bit initialises the excess ice in a way that roughly makes the surface elevation make sense - not the nicest!
        temp_ice_top = surface_heights - sp.xice_depth_bs
        j_ice_interface = np.array([np.min(np.arange(ny)[std_layer_tops < temp_ice_top[i]]) for i in range(nx)])
        xice_point = np.zeros_like(yy_temp_o)
        for i in range(right_slope):
            j = j_ice_interface[i]
            elevation_remaining = surface_heights[i]
            while elevation_remaining > 0 and j < ny:
                if elevation_remaining >= std_y_spacings[j] * sp.xice_content:
                    xice_point[i,j] = sp.xice_content
                    elevation_remaining -= std_y_spacings[j] * sp.xice_content
                    j += 1
                else: # last one gets a bit less.
                    xice_point[i,j] = elevation_remaining / std_y_spacings[j]
                    elevation_remaining = 0.0
            if elevation_remaining < 0.0 - small_value or j >= ny:
                print(f'Elevation_remaining = {elevation_remaining}, surface height target: {surface_heights[i]}, elevation so far = {np.sum(std_y_spacings*xice_point[i,:])}')
                print(f'xice deposited between layers {j_ice_interface[i]} and {j} inclusive.')
                print(f' xice in column: {xice_point[i,:]}')
                print(f'spacings * xice = {std_y_spacings*xice_point[i,:]}')
                if elevation_remaining < 0.0 - small_value:
                    raise Exception(f'Elevation_remaining < 0.0, see out log for more details.')
                if j >= ny:
                    raise Exception(f'Run out of layers for xice, see out log for more details.')

else:
    nx, ny = 20,20    
    if sp.option in ['twotilesat', 'pardry']:
        slope_y = 0.5 # transition between dry and wet 

    if sp.option == 'pardry':
        slope_x = 1.0 # transition between dry and wet
        mid = int(nx/2)
        left_slope = int(mid-(slope_x/(2*thickness)))
        right_slope = int(mid+(slope_x/(2*thickness)))
        
if sp.option != 'offset': std_x_positions = (np.arange(nx)-(nx-1)/2)*thickness
#--------------------------------------------------------------------------
# Functions
#--------------------------------------------------------------------------

def heat_con(hcon, sthu, sthf, v_sat):
    """
    Inputs must be arrays of equal shape, ignoring case of no soil points.
    Using soilhc_method=3 from heat_con_jls_mod.
    
    hcon is Dry soil thermal conductivity (W/m/K).
    """ 
    hcsat_max_allowed = 2.20
    
    thice = np.zeros(np.shape(hcon))
    thice[sthf>0] = v_sat[sthf>0] * sthf[sthf>0] / (sthu[sthf>0] + sthf[sthf>0])
        # The concentration of ice at saturation for the current mass fraction
        # of liquidcwater (m3 H2O/m3 soil).  
    thwat = v_sat - thice
        # The concentration of liquid water at saturation for the current mass
        # fraction of liquid (m3 H2O/m3 soil).
    sth = sthu + sthf
        # Fractional saturation of water (liquid+ice) at layer boundaries.
    hcsat = ( 1.0 - 0.0134 * np.log(hcon) ) / ( -0.745 - np.log(hcon) )
        # The thermal conductivity of the saturated  soil at current ratio of ice
        # to liquid water (W/m/K).
        
    hcsat[hcsat > hcsat_max_allowed] = hcsat_max_allowed
    hcsat[hcsat < 0.5]            = 0.5
    # Adjust HCSAT for frozen soil water
    hcsat = hcsat * (hcwat**thwat) * (hcice**thice) / (hcwat**v_sat)

    ke = np.zeros(np.shape(hcon))
    ke[sth > 0.1] = 1.0 + np.log10(sth[sth > 0.1]) # Kersten number

    hcons = (hcsat - hcon) * ke + hcon
    
    return hcons
        # The thermal conductivity between adjacent layers including effects of
        # water and ic (W/m/K).
        
def calc_tmax(v_sat, smcl, smclsat, bexp, sathh):
    """
    Calculate TMAX, the temperature above which all soil water is
    unfrozen.
    """
    tmax = -zerodegc*np.ones(np.shape(v_sat)) # Hang on, was this originally written for kelvin or celsius?
    index = (v_sat > tiny_0) * (smcl > small_value)
    work1 = np.ones(np.shape(v_sat)) * -999
    work1[index] = (smcl[index] / smclsat[index])**(bexp[index])
    temp_tmax = np.ones(np.shape(v_sat)) * -999
    temp_tmax[work1>small_value] = -sathh[work1>small_value] / (dpsidt * work1[work1>small_value])
    tmax = np.maximum(temp_tmax,tmax)   
    return tmax

def depth_unfrozen(sample):
    if sample.min()<0:
        return sample[sample<0][0].depth.data
    else:
        return np.nan
    
def process_input(ds, n_tasks, task_n, subtask_n, n_subtasks, temp_adjust, point_run, point_lat, point_lon, subchunk):
    ds = ds.copy().where(ds.latitude>50,drop=True) # probably already done
    ds = ds.assign_coords(land=('x',ds.x.data))
    if point_run:
        ds = ds.where((ds.latitude==point_lat)*(ds.longitude==point_lon),drop=True)
    else:
        if sp.crop in ['lon', 'latlon']:
            ds = ds.where(ds.longitude>sp.minlon,drop=True)
            ds = ds.where(ds.longitude<sp.maxlon,drop=True)
            if sp.crop == 'latlon':
                ds = ds.where(ds.latitude>sp.minlat,drop=True)
                ds = ds.where(ds.latitude<sp.maxlat,drop=True)
        elif sp.crop == 'land':
            n_landpoints = len(ds.x)
            dividers = np.linspace(0, n_landpoints, n_tasks+1).astype(int)
            task_bounds = [(dividers[i], dividers[i+1]-1) for i in range(len(dividers)-1)]
            min_land, max_land = task_bounds[task_n]
            print(f'There are {n_landpoints} landpoints, task {task_n} out of {n_tasks} has min_land = {min_land}, max_land = {max_land}')
            if subchunk:
                subdividers = np.linspace(min_land,max_land+1,n_subtasks+1).astype(int)
                subtask_bounds = [(subdividers[i], subdividers[i+1]-1) for i in range(len(subdividers)-1)]
                min_land, max_land = subtask_bounds[subtask_n]
                print(f'but, subtask {subtask_n} out of {n_subtasks} has min_land = {min_land}, max_land = {max_land}')
            ds = ds.where(ds.land>=min_land, drop=True)
            ds = ds.where(ds.land<=max_land, drop=True)

    if temp_adjust:
        raise Exception("temp_adjust broken! Where's the file folder?")
        folder = None
        hist = xa.open_dataset(folder+'gfdl-esm4_r1i1p1f1_w5e5_historical_tasAdjust_global_annual.nc')
        ssp_585 = xa.open_dataset(folder+'gfdl-esm4_r1i1p1f1_w5e5_ssp585_tasAdjust_global_annual.nc')
        hist_data= hist.tasAdjust[-5,:,:] -273.15 # MAT for 2010
        future_data = ssp_585.tasAdjust[-10:-1,:,:].mean(axis=0) -273.15 # MAT for 2090 to 2100
        temp_diff = future_data - hist_data
        temp_adjust = np.zeros(len(high['t_soil'].longitude[0,:]))
        for l in high.x.data:
            temp_adjust[l] = temp_diff.sel({'lat':high.latitude.data[0,l],'lon':high.longitude.data[0,l]}).data
        ds['t_soil'][:,0,0,:] = ds['t_soil'][:,0,0,:] + temp_adjust

    return ds

def initialise_input(task_n, n_tasks, subtask_n, n_subtasks, file_prefix, point_run, point_lat, point_lon, subchunk):
    #-----------------------------------------------------------------------------------
    # Driving data
    #-----------------------------------------------------------------------------------
    print('Loading driving data')
    if sp.duo:
        duo_folder = f'/home/users/eleanorburke/eleanorburke/noah/pfduo/GFDL-ESM4/{sp.hist_or_sspno}/'
        duo_physmon = f'isimip3b_duo_50Npfduo_gfdl-esm4_{sp.hist_or_sspno}.phys_mon'
        if sp.hist_or_sspno == 'historical':
            y_start, y_end = 1912, 2014 # 1911 and 2015 were dodgy, otherwise could have gone 1850 to 2015
        else:
            y_start, y_end = 2015, 2088
        xds_duo = xa.open_mfdataset([duo_folder + duo_physmon + f'.{year}.nc' for year in range(y_start,y_end)]) # 1850-1911, 1912-2015 (not including ends) 2015, 2088
        high = xds_duo.copy().where(xds_duo.x%2==0,drop=True)
        low = xds_duo.copy().where(xds_duo.x%2==1,drop=True)
        if sp.surf_daily:
            duo_surf_daily = 'isimip3b_duo_50Npfduo_gfdl-esm4_historical.tsoilsurf_daily'
            xds_duo_daily = xa.open_mfdataset(duo_folder + duo_surf_daily + f'.1???.nc')
            high_daily = xds_duo_daily.copy().where(xds_duo_daily.x%2==0,drop=True)
            low_daily = xds_duo_daily.copy().where(xds_duo_daily.x%2==1,drop=True)
    else:
        high_option = 'pf'
        low_option = 'pfsat'
        hist_physmon = '_gfdl-esm4_historical.phys_mon'
        hist_tsoilsurf = '_gfdl-esm4_historical.tsoilsurf_daily'
        driving_folder = '/home/users/eleanorburke/eleanorburke/noah/may2022/'
        high = xa.open_mfdataset(f'{driving_folder}isimip3b_50N{high_option}{hist_physmon}.200*.nc')
        low = xa.open_mfdataset(f'{driving_folder}isimip3b_50N{low_option}{hist_physmon}.200*.nc')
        if sp.surf_daily:
            high_daily = xa.open_mfdataset(f'{driving_folder}isimip3b_50N{high_option}{hist_tsoilsurf}.200*.nc')
            low_daily = xa.open_mfdataset(f'{driving_folder}isimip3b_50N{low_option}{hist_tsoilsurf}.200*.nc')
    
    high = process_input(high, n_tasks, task_n, subtask_n, n_subtasks, sp.temp_adjust, point_run, point_lat, point_lon, subchunk)
    low = process_input(low, n_tasks, task_n, subtask_n, n_subtasks, sp.temp_adjust, point_run, point_lat, point_lon, subchunk)
    if sp.surf_daily:
        high_daily = process_input(high_daily, n_tasks, task_n, subtask_n, n_subtasks, sp.temp_adjust, point_run, point_lat, point_lon, subchunk)
        low_daily = process_input(low_daily, n_tasks, task_n, subtask_n, n_subtasks, sp.temp_adjust, point_run, point_lat, point_lon, subchunk)
    else:
         high_daily, low_daily = None, None              

    #-----------------------------------------------------------------------------------
    # Soil properties
    #-----------------------------------------------------------------------------------
    print('Regridding global soil properties')

    props_to_copy = ['sm_sat','b','sathh','hcap','hcon']

    if sp.soil_prop_to_use == 'dynamic_dump': # Broken!
        raise Exception('Broken soil prop routine!')
        xds_duo_dump = xa.open_dataset(duo_folder + 'isimip3b_duo_dynamicsoil_gfdl-esm4_historical.dump.18620101.0.nc')
        for var in props_to_copy:
            high[var]=high['t_soil'][:,:,:,:].copy()
            for m in range(12):
                xds_duo_ilamb50[var].data[m,:,0,:] = xds_duo_dump[var].data

    elif sp.soil_prop_to_use == 'organic':
        global_org_props = xa.load_dataset('/home/users/eleanorburke/eleanorburke/organic.layeredsoil.latlon_fixed.nc')
        soil_depths = global_org_props.soil.data
        dzsoil_io = []
        new_glob_org = global_org_props.soil.data
        for i in range(len(new_glob_org)):
            if i == 0:
                dzsoil_io.append(new_glob_org[0]*2)
            else:
                dzsoil_io.append((new_glob_org[i]-new_glob_org[i-1])*2
                                - dzsoil_io[-1])
        for prop in props_to_copy:
            print(prop)
            high[prop+'_org'] = high['t_soil'][0,:,:,:].copy()
            for l in tqdm(high.x.data):
                high[prop+'_org'][:,0,l] = global_org_props.sel({'lat':high.latitude.data[0,l],'lon':high.longitude.data[0,l]})[prop].data
                                      
    elif sp.soil_prop_to_use == 'peat':
        iskoras_props = xa.load_dataset( '/home/users/nds211/full_PFsites_copier/copied/PFsites/ancillaries_v4/Iskoras/xy_soil_properties_file.nc')
        dzsoil_io=[0.05,0.08408964,0.11397535,0.14142136,0.16718508,0.19168293,
             0.21517585,0.23784142,0.25980762,0.28117066,0.30200527,
             0.32237098,0.34231625,0.36188121]+6*[0.8]
        soil_depths = [dzsoil_io[0]/2]
        for i in range(1,len(dzsoil_io)): # check this!!!!!!!!!
            soil_depths.append(soil_depths[i-1] + (dzsoil_io[i-1] + dzsoil_io[i])/2)
        for prop in props_to_copy:
            print(prop)
            high[prop+'_org'] = high['t_soil'][0,:,:,:].copy()
            high[prop+'_org'] = high[prop+'_org'].transpose('x','y','soil')
            high[prop+'_org'][:,0,:] = iskoras_props[prop][:,0,0].data
            high[prop+'_org'] = high[prop+'_org'].transpose('soil','y','x')
    
    high['soil'], low['soil'] = soil_depths, soil_depths
    high['soil'].attrs['units'], low['soil'].attrs['units'] = 'm', 'm'

    #-----------------------------------------------------------------------------------
    # Initialising
    #-----------------------------------------------------------------------------------
    print('Initialising')
    
    nl = len(high.x)
    grid = np.ones(shape=(nx,ny,nl))
    grid_to_cube = lambda xy_like: (xy_like * np.ones((nl,nx,ny))).transpose(1,2,0)
    column_to_cube = grid_to_cube # Since it works for both!
    if sp.option == 'offset':
        soil_depths = std_y_positions # grid positions, overwrites ancillary soil_depths
        dx, dy = [(x*np.ones((nl,nx,ny))).transpose(1,2,0) for x in np.meshgrid(std_x_spacings, std_y_spacings,indexing='ij')] # old way using column arrays  
        z_depth = grid_to_cube(xy_mask*soil_depths)
        grid = grid_to_cube(xy_mask)
        
        n_unsat = -1
        while np.isnan(xy_mask[-1,n_unsat+1]):
            n_unsat += 1
    else:
        dx = grid * thickness # shape of each cell
        dy = grid.copy()
        z_depth = column_to_cube(soil_depths)
        
        n_unsat = 0
        if sp.option in ['twotilesat', 'pardry']:
            slope_y = 0.5 # transition between dry and wet 
            while soil_depths[n_unsat] < slope_y:
                n_unsat += 1

    is_surface = find_surface_cells(z_depth)

    v_sat = grid.copy()
    bexp = grid.copy()
    sathh = grid.copy()
    hcap = grid.copy()
    smclsat = grid.copy()
    hcon = grid.copy()
    xice = grid.copy()

    if sp.option == 'offset':
        interp_soil_props = lambda ds : ds[:,0,:].interp(soil=-depth_bs[j,:], kwargs=dict(fill_value='extrapolate')).copy().to_numpy()
        for j in range(nx):
            # dy and dx already done. Could be quicker as many interpolations the same.
            v_sat[j,:,:] = interp_soil_props(high['sm_sat_org'])
            bexp[j,:,:] = interp_soil_props(high['b_org'])
            sathh[j,:,:] = interp_soil_props(high['sathh_org'])
            hcap[j,:,:] = interp_soil_props(high['hcap_org'])
            hcon[j,:,:] = interp_soil_props(high['hcon_org'])
        if sp.excess_ice:
            xice = xice * (xice_point * np.ones((nl,nx,ny))).transpose(1,2,0)
        else:
            xice = xice * 0.0
    else:
        dy_part = np.transpose(np.array(dzsoil_io)*np.ones((nl,ny)))
        temp_v_sat = high['sm_sat_org'][:,0,:].copy().to_numpy() # v_sat[i,n] ! Volumetric soil moisture concentration at saturation (m3 H2O/m3 soil). (sm_sat?)
        temp_bexp = high['b_org'][:,0,:].copy().to_numpy() # bexp[i,n] ! Brooks & Corey exponent.
        temp_sathh = high['sathh_org'][:,0,:].copy().to_numpy() # sathh[i,n] ! Saturated soil water pressure (m).
        temp_hcap = high['hcap_org'][:,0,:].copy().to_numpy() # Soil heat capacity (J/K/m3).
        temp_hcon = high['hcon_org'][:,0,:].copy().to_numpy() # Dry soil thermal conductivity (W/m/K).

        for j in range(nx):
            dy[j,:,:] = dy_part
            v_sat[j,:,:] = temp_v_sat
            bexp[j,:,:] = temp_bexp
            sathh[j,:,:] = temp_sathh
            hcap[j,:,:] = temp_hcap
            hcon[j,:,:] = temp_hcon

        xice = xice * 0.0
            
    if sp.option == 'offset':
        apply_mask = lambda cube_var: (xy_mask*cube_var.transpose(2,0,1)).transpose(1,2,0)
        dy      = apply_mask(dy)
        dx      = apply_mask(dx)
        v_sat   = apply_mask(v_sat)
        bexp    = apply_mask(bexp)
        sathh   = apply_mask(sathh)
        hcap    = apply_mask(hcap)
        hcon    = apply_mask(hcon)
        xice    = apply_mask(xice)

        if sp.mineral_layer:
            temp_mineral_top = surface_heights - sp.mineral_depth_bs
            j_mineral_interface = np.array([np.min(np.arange(ny)[std_layer_tops < temp_mineral_top[i]]) for i in range(nx)]) # should the 0.05 be here or before?
            sand, silt, clay = 0.01, 0.25, 0.74 #  find! maybe add organic too.
            # k_s = 10**(-2.75 - 0.64*clay + 1.26*sand) # no water flows used.
            v_sat_mineral = 0.505 - 0.037*clay - 0.142*sand
            bexp_mineral = 3.10 + 15.70*clay - 0.3*sand
            sathh_mineral = 0.01 * 10**(2.17 - 1.58*sand - 0.63*clay)
            hcap_mineral = (1-v_sat_mineral) * (clay*2.373E+6 + sand*2.133E+6 + silt*2.133E+6)
            hcon_mineral = 0.025**v_sat_mineral * 1.16025**((1-v_sat_mineral)*clay) * 1.57025**((1-v_sat_mineral)*sand) * 1.57025**((1-v_sat_mineral)*silt)
            for i in range(nx):
                v_sat[i,j_mineral_interface[i]:] = v_sat_mineral
                bexp[i,j_mineral_interface[i]:] = bexp_mineral
                sathh[i,j_mineral_interface[i]:] = sathh_mineral
                hcap[i,j_mineral_interface[i]:] = hcap_mineral
                hcon[i,j_mineral_interface[i]:] = hcon_mineral
            
    smclsat = rho_water*dx*(1-xice)*dy*v_sat # <!C22D!> # smclsat[i,n] ! The saturation moisture content of each cell (kg/m).

    # Initialise
    #------------
    tsl = grid.copy()*0.0
    smcl = smclsat.copy() # <!C22D!> # smcl[i,n] ! Soil moisture content of each cell (kg/m).
    
    if sp.option == 'offset':
        interp_vars = lambda ds, j : ds[0,:,0,:].interp(soil=-depth_bs[j,:], kwargs=dict(fill_value='extrapolate')).data
        tsl[:left_slope,:,:] = interp_vars(high.t_soil, 0) - 273.15
        tsl[right_slope+5:,:,:] = interp_vars(low.t_soil, -1) - 273.15
        if not sp.excess_ice: # hopefully still works
            for j in range(left_slope,right_slope):
                T_high = interp_vars(high.t_soil, j) - 273.15
                T_low = interp_vars(low.t_soil, j) - 273.15
                tsl[j,:,:] = T_high + ((j-left_slope)/(right_slope-left_slope))*(T_low-T_high) # linearly interpolate
        else: # add border to ensure slope doesn't immediately thaw
            for j in range(left_slope,right_slope):
                T_high = (interp_vars(high.t_soil, j) - 273.15).compute()
                tsl[j,:,:] = T_high
            for j in range(right_slope,right_slope+5):
                T_high = (interp_vars(high.t_soil, j) - 273.15).compute()
                T_low = (interp_vars(low.t_soil, j) - 273.15).compute()
                tsl[j,:,:] = T_high + ((j-right_slope)/5)*(T_low-T_high) # linearly interpolate
        if sp.grid_option == 'doule_sided':
            tsl[:int(nx/2),:,:] = tsl[int(nx/2):,:,:] # check this!
        
        interp_vars_unsat = lambda ds, j : ds[0,:,0,:].interp(soil=-depth_bs[j,:n_unsat], kwargs=dict(fill_value='extrapolate')).data
        if sp.grid_option == 'double_sided': # This needs looking at !!!!
            smcl[:,:n_unsat,:] = interp_vars_unsat(high.sthu + high.sthf, int(nx/2)) * smclsat[:,:n_unsat,:]
        else:
            smcl[:,:n_unsat,:] = interp_vars_unsat(high.sthu + high.sthf, 0) * smclsat[:,:n_unsat,:]
        # nans in smclsat hould take care of surface? -> check! have deleted : right_slope
        smcl[:,n_unsat:,:] = smclsat[:,n_unsat:,:]
        T_high_surface, T_low_surface = (high.t_soil[0,0,0,:].data- 273.15).compute(), (low.t_soil[0,0,0,:].data- 273.15).compute() # initially not shifted
        T_surf_gradient_grid = (T_low_surface-T_high_surface)*np.ones((nx,ny,nl))
        T_surf_high_grid = T_high_surface*np.ones((nx,ny,nl))

    else:
        tsl[:int(nx/2),:,:] = high.t_soil.data[0,:,0,:] - 273.15 # Soil cell temperatures need to be in Celsius!
        tsl[int(nx/2):,:,:] = low.t_soil.data[0,:,0,:] - 273.15 # tsl[i,n] ! Soil cell temperatures (Celsius)

        if sp.option == 'twotilesat':
            smcl[0,:n_unsat,:] = (high.sthu.data[0,:n_unsat,0,:] + high.sthf.data[0,:n_unsat,0,:])*smclsat[0,:n_unsat,:]
            tsl[0,:,:] = high.t_soil.data[0,:,0,:] - 273.15 # Soil cell temperatures need to be in Celsius!
            tsl[1,:,:] = low.t_soil.data[0,:,0,:] - 273.15 # tsl[i,n] ! Soil cell temperatures (Celsius)
        else:
            tsl[:int(nx/2),:,:] = high.t_soil.data[0,:,0,:] - 273.15 # Soil cell temperatures need to be in Celsius!
            tsl[int(nx/2):,:,:] = low.t_soil.data[0,:,0,:] - 273.15 # tsl[i,n] ! Soil cell temperatures (Celsius)

            if sp.option == 'high_and_dry':
                smcl[:int(nx/2),:,:]=(high.sthu.data[0,:,0,:]+high.sthf.data[0,:,0,:])*smclsat[:int(nx/2),:,:]

            elif sp.option == 'pardry':
                T_high = high.t_soil.data[0,:,0,:] - 273.15 # Soil cell temperatures need to be in Celsius!
                T_low = low.t_soil.data[0,:,0,:] - 273.15
                smcl[:right_slope,:n_unsat,:] =(high.sthu.data[0,:n_unsat,0,:] + high.sthf.data[0,:n_unsat,0,:]) * smclsat[:right_slope,:n_unsat,:]
                for j in range(left_slope,right_slope):
                    nint = int(n_unsat*(1-((j-left_slope)/(right_slope-left_slope))))
                    smcl[j,nint:,:]=smclsat[j,nint:,:]
                    tsl[j,:,:] = T_high + ((j-left_slope)/(right_slope-left_slope))*(T_low-T_high) # linearly interpolate

    tmax = calc_tmax(v_sat, smcl, smclsat, bexp, sathh)

    # Initialise soil moistures
    # -----------------------------
    smclu = grid.copy() * 0.0
    smclf = grid.copy() * 0.0
    smclu[tsl >= tmax] = smcl[tsl >= tmax] # <!C22D!> # Unfrozen moisture content of each cell (kg/m).
    smclf[tsl >= tmax] = 0.0 # <!C22D!> # Frozen moisture content of each cell (kg/m).

    smclu[tsl<tmax] = smclsat[tsl<tmax] * (-dpsidt * tsl[tsl<tmax] / sathh[tsl<tmax])**(-1.0 / bexp[tsl<tmax])
    smclf[tsl<tmax] = smcl[tsl<tmax] - smclu[tsl<tmax] # Hmm, what if smcl < smclu?!

    sthu  = smclu / smclsat # Fractional saturation of unfrozen water.
    sthf  = smclf / smclsat # Fractional saturation of frozen water.

    #-----------------------------------------------------------------------------------
    # Save setup dump
    #-----------------------------------------------------------------------------------
    print(f'Saving setup dump')
    
    soil_prop_dump = xa.Dataset(
        data_vars = {
            'dx':(['x','depth','land'],dx),
            'dy':(['x','depth','land'],dy),
            'v_sat':(['x','depth','land'],v_sat),
            'bexp':(['x','depth','land'],bexp),
            'sathh':(['x','depth','land'],sathh),
            'hcap':(['x','depth','land'],hcap),
            'hcon':(['x','depth','land'],hcon),
            'smclsat':(['x','depth','land'],smclsat),
            'xice':(['x','depth','land'],xice)},
         coords={
            'depth': soil_depths, 
            'x': std_x_positions,
            'land':high['land'].data,
            'latitude':('land',high['latitude'][0,:].data),
            'longitude':('land',high['longitude'][0,:].data)},
         attrs={
            'depth':'Depth (m)',
            'x':'x (m)',
            'option':sp.option,
            'duo':str(sp.duo),
            'surf_daily':str(sp.surf_daily),
            'soil_prop_to_use':sp.soil_prop_to_use,
            'sp':str(vars(sp))})
    for var, long_name, unit in [('dx','Cell height','m'),
                                 ('dy','Cell width', 'm'),
                                 ('v_sat','Volumetric soil moisture concentration at saturation', 'm3 H2O/m3 soil'),
                                 ('bexp','Brooks & Corey exponent', 'm'),
                                 ('sathh','Saturated soil water pressure', 'm'),
                                 ('hcap','Soil heat capacity', 'J/K/m3'),
                                 ('hcon','Dry soil thermal conductivity', 'W/m/K'),
                                 ('smclsat','Cell saturation moisture content', 'kg/m'),
                                 ('xice','Volumetric fraction of non-pore ice', 'm3/m3')]:
        soil_prop_dump[var].attrs['long_name'] = long_name
        soil_prop_dump[var].attrs['unit'] = unit
    
    pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in soil_prop_dump.data_vars}
    soil_prop_dump.to_netcdf( output_path + file_prefix+f'.setup_dump.nc', encoding=encoding)
    
    return nl, dx, dy, soil_depths, high, low, high_daily, low_daily, v_sat, bexp, sathh, hcap, hcon, smclsat, tsl, smcl, tmax, smclu, smclf, sthu, sthf, n_unsat, xice, z_depth, is_surface, T_surf_gradient_grid, T_surf_high_grid
                        
def initialise_outputs(output_profiles, end_datetime, start_datetime, soil_depths,
                      nx, ny, nl, thickness, high):
    """ Returns outputs and total_timedeltas. """
    outputs = output_profiles.copy() # in parm mod

    total_timedeltas = {
        'days': (end_datetime - start_datetime).days,
        'months': 12 * (end_datetime.year - start_datetime.year) \
            + (end_datetime.month - start_datetime.month),
        'years': end_datetime.year - start_datetime.year}

    for output in outputs:
        # ignoring multiple runs, timestamps now end of day / month - are years ok?
        n_timeout, timestamps = {'daily':(total_timedeltas['days'],
                                          np.arange(np.datetime64(start_datetime),np.datetime64(end_datetime), 
                                                    np.timedelta64(1, 'D')).astype('datetime64[D]') + pd.Timedelta(days=1) ),
                             'monthly':(total_timedeltas['months'],
                                        (np.arange(np.datetime64(start_datetime,'M'), np.datetime64(end_datetime, 'M'),
                                                  np.timedelta64(1, 'M'))+np.timedelta64(1, 'M')).astype('datetime64[D]')),
                             'yearly':(total_timedeltas['years'],
                                       np.arange(np.datetime64(start_datetime,'Y'), np.datetime64(end_datetime, 'Y'),
                                                 np.timedelta64(1, 'Y')).astype('datetime64[D]'))}[output['freq']]
        # Not saving initial conditions, saving at the end of a day / month / year etc.
        # Could even output as grid...
        
        if output['vars'] == 'ends':
            nx_output = 2
            output_x_coords = np.array([-nx*thickness/2,nx*thickness/2])
        else: 
            nx_output = nx
            output_x_coords = std_x_positions
            
        possible_vars = {'t_soil':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                         'sthu':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'sthf':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'sthuf':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'lat_thaw':(['time','land'],np.full((n_timeout,nl),np.nan)),
                        'lat_thaw_layers':(['time','depth','land'],np.full((n_timeout,ny,nl),np.nan)),
                        'depth_unfrozen':(['time','x','land'],np.full((n_timeout,2,nl),np.nan)),
                        'hcapt': (['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'k_y_grid':(['time','x','layers','land'],np.full((n_timeout,nx_output,ny-1,nl),np.nan)),
                        'remained_frozen':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'xice':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'xice_total':(['time','x','land'],np.full((n_timeout,nx_output,nl),np.nan)),
                        'z_depth_surface':(['time','x','land'],np.full((n_timeout,nx_output,nl),np.nan)),
                        'dy':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'tmax':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'smcl':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'smclsat':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'v_sat':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'bexp':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'hcap':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'hcon':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'sathh':(['time','x','depth','land'],np.full((n_timeout,nx_output,ny,nl),np.nan)),
                        'x_of_y_threshold_lower':(['time','land'],np.full((n_timeout,nl),np.nan)),
                        'x_of_y_threshold_upper':(['time','land'],np.full((n_timeout,nl),np.nan)),
                        'slope_gradient':(['time','land'],np.full((n_timeout,nl),np.nan)),
                        'x_offset':(['time','land'],np.full((n_timeout,nl),np.nan))}
        
        var_infos = {'t_soil' : ['Temperature','Celsius'],
                     'sthf' : ['Frozen soil moisture', 'Fractional saturation'],
                     'sthu' : ['Unfrozen soil moisture', 'Fractional saturation'],
                     'sthuf': ['Total soil moisture', 'Fractional saturation'],
                     'lat_thaw' : ['Lateral thaw', 'J/m'],
                     'lat_thaw_layers' : ['Lateral thaw', 'J/m'],
                     'depth_unfrozen' : ['Depth unfrozen','m'],
                     'hcapt' : ['Total (soil+water) volumetric heat capacity', 'J/m3/K'],
                     'k_y_grid' : ['Total (soil+water) thermal conductivity in y', 'W/m/K'],
                     'remained_frozen' : ['Soil layer remained frozen', 'Y/N'],
                     'xice' : ['Volumetric fraction of non-pore ice', 'm3/m3'],
                     'xice_total' : ['Non-pore ice soil column total', 'm3/m2'],
                     'z_depth_surface' : ['Surface elevation', 'm'],
                     'dy' : ['Cell thickness', 'm'],
                     'tmax' : ['Maximum temperature of frozen pore ice', 'Celsius'],
                     'smcl' : ['Cell moisture content', 'kg/m'],
                     'smclsat' : ['Cell saturation moisture content', 'kg/m'],
                     'v_sat' : ['Volumetric soil moisture concentration at saturation', 'm3 H2O/m3 soil'],
                     'bexp' : ['Brooks & Corey exponent', 'm'],
                     'hcap' : ['Soil heat capacity', 'J/K/m3'],
                     'hcon' : ['Dry soil thermal conductivity', 'W/m/K'],
                     'sathh' : ['Saturated soil water pressure', 'm'],
                     'x_of_y_threshold_lower' : [f'x position of lower ({sp.y_threshold_lower}m) surface elevation threshold', 'm'],
                     'x_of_y_threshold_upper' : [f'x position of upper ({sp.y_threshold_upper}m) surface elevation threshold', 'm'],
                    'slope_gradient' : [f'Gradient between surface elevation thresholds at {sp.y_threshold_lower} and {sp.y_threshold_upper}m', 'm/m'],
                     'x_offset' : ['shifting domian x offset', 'm']}
        
        output['data']=xa.Dataset(
                data_vars = {},
                 coords={
                    'time':timestamps,
                    'depth': soil_depths, 
                    'x': output_x_coords,
                    'land':high['land'].data,
                    'latitude':('land',high['latitude'][0,:].data),
                    'longitude':('land',high['longitude'][0,:].data)},
                 attrs={
                    'time': 'Date',
                    'depth':'Depth (m)',
                    'x':'x (m)'})
        
        var_profiles = {'full' : ['t_soil', 'sthf', 'sthuf'], 
                        'reduced_year' : ['lat_thaw', 'lat_thaw_layers','remained_frozen'],
                        'ends' : ['t_soil', 'sthf', 'sthuf', 'depth_unfrozen'],
                        'debug' : ['t_soil', 'sthf', 'sthuf', 'hcapt', 'k_y_grid'],
                        'debug_xice' : ['t_soil', 'sthf', 'sthuf', 'hcapt', 'k_y_grid','tmax','smcl','smclsat','v_sat','bexp', 'xice'] ,
                        'min_xice':[]} # 'lat_thaw', 'lat_thaw_layers','remained_frozen' were in full
        
        if sp.excess_ice:
            for profile in var_profiles:
                if profile == 'full':
                    var_profiles[profile].append('xice')
                var_profiles[profile].append('xice_total')
                if profile != 'ends':
                    var_profiles[profile].append('z_depth_surface')
                    var_profiles[profile].append('x_of_y_threshold_lower')
                    var_profiles[profile].append('x_of_y_threshold_upper')
                    var_profiles[profile].append('slope_gradient')
                    if sp.shifting_domain:
                        var_profiles[profile].append('x_offset')
        
        for var in var_profiles[output['vars']]:
            output['data'][var] = possible_vars[var]
            output['data'][var].attrs['long_name']=var_infos[var][0]
            output['data'][var].attrs['unit']=var_infos[var][1]
            
        comp = dict(zlib=True, complevel=5, _FillValue=None)
        output['encoding'] = {var: comp for var in output['data'].data_vars}
        
    return outputs, total_timedeltas

def update_driving(previous_time, current_time, tsl, smcl, high, 
                   low, smclsat, tmax, high_daily, low_daily, 
                   n_unsat, v_sat, bexp, sathh, z_depth, is_surface, T_surf_gradient_grid,T_surf_high_grid, nl):
    """
    Possibly updates the driving data (returns tsl, smcl and tmax)
    
    Caution: option, nx, right_slope, left_slope not included in args.
    
    """
    #month = current_time.month # this outputs the month as an int - not what we want.
    
    if sp.surf_daily:
         if current_time.day != previous_time.day:
            current_date = np.datetime64(current_time,'D')
            T_high = high_daily.t_soil_surf.loc[current_date,0,:] - 273.15
            T_low = low_daily.t_soil_surf.loc[current_date,0,:] - 273.15
            if sp.option == 'offset':
                tsl[:left_slope,0,:] = T_high
                tsl[right_slope:,n_unsat+1,:] = T_low
                for j in range(left_slope,right_slope):
                    tsl[j,surf_y_index[j],:] = T_high + ((j-left_slope)/(right_slope-left_slope))*(T_low-T_high)
            else:
                tsl[:int(nx/2),0,:] = T_high
                tsl[int(nx/2):,0,:] = T_low
                if sp.option == 'pardry':
                    for j in range(left_slope,right_slope):
                        tsl[j,0,:] = T_high + ((j-left_slope)/(right_slope-left_slope))*(T_low-T_high) # linearly interpolate surface temp
#             elif option == 'normale': # from when this was also daily
#                 smcl[:int(nx/2),:,:]=(high.sthu.loc[current_date,:,0,:].data+high.sthf.loc[current_date,:,0,:].data)*smclsat[:int(nx/2),:,:]
#                 smcl[int(nx/2):,:,:]=(low.sthu.loc[current_date,:,0,:].data+low.sthf.loc[current_date,:,0,:].data)*smclsat[int(nx/2):,:,:]
#                 tmax = calc_tmax(v_sat, smcl, smclsat, bexp, sathh) # as updating smcl

    if current_time.month != previous_time.month: # iterate monthly driving data
        driving_month = np.datetime64(current_time,'M') + np.timedelta64(1, 'M') # as timestamp at end, and current time is still the old month. Outputs date rounded to month.
        if sp.option == 'offset':
            interp_vars_unsat = lambda ds, j : ds.loc[driving_month,:,0,:].interp(soil=-depth_bs[j,:n_unsat], kwargs=dict(fill_value='extrapolate')).data
            if not sp.excess_ice:
                smcl[:,:n_unsat,:] = interp_vars_unsat(high.sthu + high.sthf, 0) * smclsat[:,:n_unsat,:]
            else:
                if sp.grid_option == 'double_sided': # This needs looking at !!!!
                    smcl[:,:n_unsat,:] = interp_vars_unsat(high.sthu + high.sthf, int(nx/2)) * smclsat[:,:n_unsat,:]
                else:
                    smcl[:,:n_unsat,:] = interp_vars_unsat(high.sthu + high.sthf, 0) * smclsat[:,:n_unsat,:]
                # Do I need to mask this? Seem not to have previously!
                # Also, just using layer depths rather than accounting for changes in surface depth.
            
        if sp.option == 'twotilesat':
            tsl[0,0,:] = high.t_soil.loc[driving_month,:,0,:][0,:].data - 273.15
            tsl[1,0,:] = low.t_soil.loc[driving_month,:,0,:][0,:].data - 273.15
            smcl[0,:n_unsat,:]=(high.sthu.loc[driving_month,:,0,:][:n_unsat,:].data+high.sthf.loc[driving_month,:,0,:][:n_unsat,:].data)*smclsat[0,:n_unsat,:]
        else:
            T_high = high.t_soil.loc[driving_month,:,0,:][0,:].data - 273.15 + sp.T_shift
            T_low = low.t_soil.loc[driving_month,:,0,:][0,:].data - 273.15 + sp.T_shift
            T_high_surface = T_high.compute()
            T_low_surface = T_low.compute()
            T_surf_gradient_grid = (T_low_surface-T_high_surface)*np.ones((nx,ny,nl))
            T_surf_high_grid = T_high_surface*np.ones((nx,ny,nl))

            if sp.option == 'offset' and not sp.excess_ice:
                tsl[:left_slope,0,:] = T_high
                tsl[right_slope:,n_unsat+1,:] = T_low
                for j in range(left_slope,right_slope):
                    tsl[j,surf_y_index[j],:] = T_high + ((j-left_slope)/(right_slope-left_slope))*(T_low-T_high)
            elif sp.option == 'offset' and sp.excess_ice:
                tsl[is_surface] = T_surf_gradient_grid[is_surface]*(sp.high_bit_elevation - std_y_spacings[0]/2 - z_depth[is_surface]) / sp.high_bit_elevation + T_surf_high_grid[is_surface] # Assuming surface layer depths the same in driving and here
            else:
                if not sp.surf_daily:
                    tsl[:int(nx/2),0,:] = T_high
                    tsl[int(nx/2):,0,:] = T_low
                if not 'no_heat_flux_edges':
                    tsl[0,:,:] = high.t_soil.loc[driving_month,:,0,:].data - 273.15
                    tsl[-1,:,:] = low.t_soil.loc[driving_month,:,0,:].data - 273.15
                    tsl[:int(nx/2),-1,:] = high.t_soil.loc[driving_month,:,0,:][-1,:].data - 273.15
                    tsl[int(nx/2):,-1,:] = low.t_soil.loc[driving_month,:,0,:][-1,:].data - 273.15
                if sp.option == 'high_and_dry':
                    smcl[:int(nx/2),:,:]=(high.sthu.loc[driving_month,:,0,:].data+high.sthf.loc[driving_month,:,0,:].data)*smclsat[:int(nx/2),:,:]
                elif sp.option == 'pardry':
                    smcl[:right_slope,:n_unsat,:] =(high.sthu.loc[driving_month,:,0,:][:n_unsat,:].data + high.sthf.loc[driving_month,:,0,:][:n_unsat,:].data)*smclsat[:right_slope,:n_unsat,:] # why not index by depth instead of finding n_unsat?
                    for j in range(left_slope,right_slope):
                        nint = int(n_unsat*(1-((j-left_slope)/(right_slope-left_slope))))
                        smcl[j,nint:,:]=smclsat[j,nint:,:]
                        if not sp.surf_daily:
                            tsl[j,0,:] = T_high + ((j-left_slope)/(right_slope-left_slope))*(T_low-T_high) # linearly interpolate surface temp
                elif sp.option == 'normale': # to test
                    if not sp.surf_daily:
                        T_high = high.t_soil.loc[driving_month,:,0,:][0,:].data - 273.15 # Soil cell temperatures need to be in Celsius!
                        T_low = low.t_soil.loc[driving_month,:,0,:][0,:].data - 273.15
                    smcl[:int(nx/2),:,:]=(high.sthu.loc[driving_month,:,0,:].data+high.sthf.loc[driving_month,:,0,:].data)*smclsat[:int(nx/2),:,:]
                    smcl[int(nx/2):,:,:]=(low.sthu.loc[driving_month,:,0,:].data+low.sthf.loc[driving_month,:,0,:].data)*smclsat[int(nx/2):,:,:]
        tmax = calc_tmax(v_sat, smcl, smclsat, bexp, sathh) # only updating smcl monthly
                        
    return tsl, smcl, tmax, T_surf_gradient_grid, T_surf_high_grid

def interp_soilvar_faster(soilvar_eff, z_eff, zsoil_extend):
    """
    Takes the effective vertical soil variable soilvar_eff[dim_cslayer, dim_soilvar] 
    on the effective zsoil depths (m) z_eff[dim_cslayer], 
    and returns it reinterpolated as soilvar[dim_cslayer, dim_soilvar]
    onto the original zsoil depths (with extended base layer) zsoil_extend[dim_cslayer].
    
    dim_cslayer is the number of vertical soil carbon layers and
    dim_soilvar is the number of pools in soil var.

    Note: here, we are using the dim_soilvar dimension to input multiple soil vars at once.

    Using Chadburn et al. (2022) - https://doi.org/10.5194/gmd-15-1633-2022
    """
    dim_soilvar, dim_cslayer = np.shape(soilvar_eff)
    # Interpolate the carbon densities to the correct layer points
    depth_term_plus = np.zeros(dim_cslayer)
    depth_term_minus = np.zeros(dim_cslayer) 
    depth_term_plus[:-1] = (z_eff[:-1] - zsoil_extend[:-1]) / (z_eff[:-1] - z_eff[1:]) # see other plus, is this consistent?
    depth_term_minus[1:] = (z_eff[1:] - zsoil_extend[1:]) / (z_eff[1:] - z_eff[:-1]) # see other minus, is this consistent?

    # sort out multipliers
    mult_75 = np.full(dim_cslayer, 0.75)
    mult_75[0] = 1.0
    mult_75[-1] = 0.0   
    mult_75_plus = np.zeros((dim_soilvar, dim_cslayer)) # not original, neccessary?
    mult_75_minus = np.zeros((dim_soilvar, dim_cslayer)) # not original, neccessary?
    mult_75_plus[:,:-1] = mult_75[:-1] * (soilvar_eff[:,1:] - soilvar_eff[:,:-1])
    mult_75_minus[:,1:] = mult_75[1:] * (soilvar_eff[:,:-1] - soilvar_eff[:,1:])

    mult_25 = np.full(dim_cslayer, 0.25)
    mult_25[0] = 0.0
    mult_25[-1] = 1.0
    mult_25_plus = np.zeros((dim_soilvar, dim_cslayer)) # not original, neccessary?
    mult_25_minus = np.zeros((dim_soilvar, dim_cslayer)) # not original, neccessary?
    mult_25_plus[:,:-1] = mult_25[:-1] * (soilvar_eff[:,1:] - soilvar_eff[:,:-1])
    mult_25_minus[:,1:] = mult_25[1:] * (soilvar_eff[:,:-1] - soilvar_eff[:,1:])

    soilvar = np.zeros((dim_soilvar, dim_cslayer))

    lower = z_eff < zsoil_extend
    lower[[0,-1]] = True

    soilvar[:,lower] = (soilvar_eff[:,lower] +
                          depth_term_plus[lower] * mult_75_plus[:,lower] +
                          depth_term_minus[lower] * mult_25_minus[:,lower])  

    soilvar[:,~lower] = (soilvar_eff[:,~lower] +
                          depth_term_minus[~lower] * mult_75_minus[:,~lower] +
                          depth_term_plus[~lower] * mult_25_plus[:,~lower])
    return soilvar

def interp_soilvar(soilvar_eff, z_eff, zsoil_extend):
    """
    Takes the effective vertical soil variable soilvar_eff[dim_cslayer, dim_soilvar] 
    on the effective zsoil depths (m) z_eff[dim_cslayer], 
    and returns it reinterpolated as soilvar[dim_cslayer, dim_soilvar]
    onto the original zsoil depths (with extended base layer) zsoil_extend[dim_cslayer].
    
    dim_cslayer is the number of vertical soil carbon layers and
    dim_soilvar is the number of pools in soil var.
    """
    dim_soilvar, dim_cslayer = np.shape(soilvar_eff)
    # Interpolate the carbon densities to the correct layer points
    depth_term_plus = np.zeros(dim_cslayer)
    depth_term_minus = np.zeros(dim_cslayer) 
    depth_term_plus[:-1] = (z_eff[:-1] - zsoil_extend[:-1]) / (z_eff[:-1] - z_eff[1:]) # see other plus, is this consistent?
    depth_term_minus[1:] = (z_eff[1:] - zsoil_extend[1:]) / (z_eff[1:] - z_eff[:-1]) # see other minus, is this consistent?

    # sort out multipliers
    mult_75 = np.full(dim_cslayer, 0.75)
    mult_75[0] = 1.0
    mult_75[-1] = 0.0   
    mult_75_plus = np.zeros((dim_soilvar, dim_cslayer)) # not original, neccessary?
    mult_75_minus = np.zeros((dim_soilvar, dim_cslayer)) # not original, neccessary?
    mult_75_plus[:,:-1] = mult_75[:-1] * (soilvar_eff[:,1:] - soilvar_eff[:,:-1])
    mult_75_minus[:,1:] = mult_75[1:] * (soilvar_eff[:,:-1] - soilvar_eff[:,1:])

    mult_25 = np.full(dim_cslayer, 0.25)
    mult_25[0] = 0.0
    mult_25[-1] = 1.0
    mult_25_plus = np.zeros((dim_soilvar, dim_cslayer)) # not original, neccessary?
    mult_25_minus = np.zeros((dim_soilvar, dim_cslayer)) # not original, neccessary?
    mult_25_plus[:,:-1] = mult_25[:-1] * (soilvar_eff[:,1:] - soilvar_eff[:,:-1])
    mult_25_minus[:,1:] = mult_25[1:] * (soilvar_eff[:,:-1] - soilvar_eff[:,1:])

    soilvar = np.zeros((dim_soilvar, dim_cslayer))
    for n in range(dim_cslayer):
        if ((n == 0 or n == dim_cslayer) or z_eff[n] < zsoil_extend[n]): # surely n will never == dim_sclayer? fortran to python issue?
            soilvar[:,n] = (soilvar_eff[:,n] +
                          depth_term_plus[n] * mult_75_plus[:,n] +
                          depth_term_minus[n] * mult_25_minus[:,n])  
        else: # z_eff >= zsoil_extend
            soilvar[:,n] = (soilvar_eff[:,n] +
                          depth_term_minus[n] * mult_75_minus[:,n] +
                          depth_term_plus[n] * mult_25_plus[:,n])
    return soilvar

def interp_soilvar_1st_order(soilvar_eff, z_eff, zsoil_extend):
    """
    A first order only version of interp_soilvar 
    - avoids all instabilities but increased smearing.

    Seems to break things.

    By default this method assumes the up direction for z is positive,
    and that input arrays are in order from the bottom to the top. 
    If they are in order from the top to the bottom, z up still positive, 
    set topdown = True.
    """
    dim_soilvar, dim_cslayer = np.shape(soilvar_eff)

    topdown = True # this could be automatic?
    if topdown: # Arrays are in order from the top to the bottom of the soil column
        soilvar_eff = soilvar_eff[:,::-1]
        z_eff = z_eff[::-1]
        zsoil_extend = zsoil_extend[::-1]

    dz = zsoil_extend - z_eff
    new_soilvar = np.zeros((dim_soilvar, dim_cslayer)) # could actually do one that can cope with nans you know.
    
    z_dz = dz == 0
    new_soilvar[:,z_dz] = soilvar_eff[:,z_dz]
    
    if dz[-1] > 0:
        new_soilvar[:,-1] = soilvar_eff[:,-1]
    p_dz = dz[:-1] > 0
    new_soilvar[:,:-1][:,p_dz] = soilvar_eff[:,:-1][:,p_dz] + dz[:-1][p_dz] * (
                        (soilvar_eff[:,1:][:,p_dz] - soilvar_eff[:,:-1][:,p_dz]) / (z_eff[1:][p_dz] - z_eff[:-1][p_dz]))
    
    if dz[0] < 0: # <
        new_soilvar[:,0] = soilvar_eff[:,0]
    n_dz = dz[1:] < 0 # <
    new_soilvar[:,1:][:,n_dz] = soilvar_eff[:,1:][:,n_dz] + dz[1:][n_dz] * (
                        (soilvar_eff[:,:-1][:,n_dz] - soilvar_eff[:,1:][:,n_dz]) / (z_eff[:-1][n_dz] - z_eff[1:][n_dz]))

    if topdown:
        soilvar_to_return = new_soilvar[:,::-1]
    else:
        soilvar_to_return = new_soilvar

    return soilvar_to_return

def find_surface_cells(z_depth): 
    """
    Surface points are in contact with the air so are currently driven by surface temperatures. True if no layers above, or if there is a blank or partially filled cell to either side.

    Uses z_depth[nx, ny, nl] to check for nan cells, returns is_surface[nx, ny, nl].
    """
    # Could change this so only checking changed points & use top_layer_index?

    is_surface = np.full(np.shape(z_depth),False)

    temp_indices = np.indices(z_depth.shape)
    y_indices = temp_indices[1].astype(float)
    y_indices[np.isnan(z_depth)] = np.nan
    y_indices = np.nanmin(y_indices, axis=1)
    y_indices = y_indices.astype(int)
    is_surface[temp_indices[0][:,0,:],y_indices,temp_indices[2][:,0,:]] = True # Check above

    is_land = ~np.isnan(z_depth)
    is_surface[:-1,:,:] += is_land[:-1,:,:] * (~is_land[1:,:,:]) # Check to right
    is_surface[1:,:,:] += is_land[1:,:,:] * (~is_land[:-1,:,:]) # Check to left

    return is_surface

# New for multiarray

def thickness_to_positions_xl_y(thicknesses = None, center = False):
    """ 
    For converting thicknesses[xl,y] of layers to midpoint positions[xl,y].
    
    """
    positions = thicknesses/2
    positions[:,1:] += np.cumsum(thicknesses[:,:-1], axis = 1)
    if center:
        positions = positions - np.sum(thicknesses, axis = 1)/2
    return positions


def interp_soilvar_multiarray(soilvar_eff, z_eff, zsoil_extend, tlis):
    """
    Takes the effective vertical soil variable soilvar_eff[dim_cslayer, dim_xl_list, dim_soilvar] 
    on the effective zsoil depths (m) z_eff[dim_xl_list, dim_cslayer], 
    and returns it reinterpolated as soilvar[dim_cslayer, dim_xl_list, dim_soilvar]
    onto the original zsoil depths (with extended base layer) zsoil_extend[dim_xl_list, dim_cslayer].
    
    dim_cslayer is the number of vertical soil carbon layers and
    dim_soilvar is the number of pools in soil var.

    Note: here, we are using the dim_soilvar dimension to input multiple soil vars at once.

    Using Chadburn et al. (2022) - https://doi.org/10.5194/gmd-15-1633-2022
    """
    dim_soilvar, dim_xl_list, dim_cslayer = np.shape(soilvar_eff)
    # Interpolate the carbon densities to the correct layer points
    depth_term_plus = np.zeros((dim_xl_list, dim_cslayer))
    depth_term_minus = np.zeros((dim_xl_list, dim_cslayer))
    depth_term_plus[:,:-1] = (z_eff[:,:-1] - zsoil_extend[:,:-1]) / (z_eff[:,:-1] - z_eff[:,1:]) # see other plus, is this consistent?
    depth_term_minus[:,1:] = (z_eff[:,1:] - zsoil_extend[:,1:]) / (z_eff[:,1:] - z_eff[:,:-1]) # see other minus, is this consistent?

    # sort out multipliers
    mult_75 = np.full((dim_xl_list, dim_cslayer), 0.75)
    mult_75[:,0] = 1.0 # hopefully not needed now
    mult_75[tlis[0],tlis[1]] = 1.0
    mult_75[:,-1] = 0.0   
    
    mult_75_plus = np.zeros((dim_soilvar, dim_xl_list, dim_cslayer)) # not original, neccessary?
    mult_75_minus = np.zeros((dim_soilvar,dim_xl_list, dim_cslayer)) # not original, neccessary?
    mult_75_plus[:,:,:-1] = mult_75[:,:-1] * (soilvar_eff[:,:,1:] - soilvar_eff[:,:,:-1])
    mult_75_minus[:,:,1:] = mult_75[:,1:] * (soilvar_eff[:,:,:-1] - soilvar_eff[:,:,1:])

    mult_25 = np.full((dim_xl_list, dim_cslayer), 0.25)
    mult_25[:,0] = 0.0
    mult_25[tlis[0],tlis[1]] = 0.0
    mult_25[:,-1] = 1.0

    mult_25_plus = np.zeros((dim_soilvar, dim_xl_list, dim_cslayer)) # not original, neccessary?
    mult_25_minus = np.zeros((dim_soilvar, dim_xl_list, dim_cslayer)) # not original, neccessary?
    mult_25_plus[:,:,:-1] = mult_25[:,:-1] * (soilvar_eff[:,:,1:] - soilvar_eff[:,:,:-1])
    mult_25_minus[:,:,1:] = mult_25[:,1:] * (soilvar_eff[:,:,:-1] - soilvar_eff[:,:,1:])

    soilvar = np.zeros((dim_soilvar, dim_xl_list, dim_cslayer))

    xl_cs_indices = np.indices((dim_xl_list, dim_cslayer))
    lower = z_eff < zsoil_extend
    lower[:,[0,-1]] = True
    lower[tlis[0],tlis[1]] = True

    lower_xl = xl_cs_indices[0][lower]
    lower_cs = xl_cs_indices[1][lower]
    soilvar[:,lower] = (soilvar_eff[:,lower_xl,lower_cs] +
                          depth_term_plus[lower_xl,lower_cs] * mult_75_plus[:,lower_xl,lower_cs] +
                          depth_term_minus[lower_xl,lower_cs] * mult_25_minus[:,lower_xl,lower_cs])  
    
    lower_xl = xl_cs_indices[0][~lower]
    lower_cs = xl_cs_indices[1][~lower]
    soilvar[:,~lower] = (soilvar_eff[:,lower_xl,lower_cs] +
                          depth_term_minus[lower_xl,lower_cs] * mult_75_minus[:,lower_xl,lower_cs] +
                          depth_term_plus[lower_xl,lower_cs] * mult_25_plus[:,lower_xl,lower_cs])

    soilvar[:,np.isnan(zsoil_extend)] = np.nan # Probably not necessary right?
    return soilvar

def excess_ice_calc_multiarray(xice, thawing_xice, dxice_tot, dy, z_depth,
                    nx, ny, nl, no_heat_flux_edges, dx,
                    std_y_positions, std_y_spacings, std_layer_tops,
                    tsl, tmax, smcl, smclu, smclf, smclsat, sthu, sthf, v_sat, bexp, sathh, hcap, hcon):
    """
    Might not work yet, and seems not to be faster!
    
    """
    if no_heat_flux_edges:
        thawing_xice[[0,nx-1],:,:] = False
    x_indices = np.indices((nx,nl))[0][thawing_xice.max(axis=1)]
    l_indices = np.indices((nx,nl))[1][thawing_xice.max(axis=1)]
    
    eff_y_spacings = dy[x_indices,:,l_indices] * (1 + dxice_tot[x_indices,:,l_indices])
    z_eff = tot_depth + thickness_to_positions_xl_y(eff_y_spacings[:,::-1])[:,::-1] # why are we reversing twice here?
    ztop = tot_depth + np.nansum(eff_y_spacings,axis = 1) # have we defined tot_depth?
    nxl = eff_y_spacings.shape[0]
    
    # Here
    temp_ycol_to_xl_y = lambda ycol_array: (ycol_array * np.ones(eff_y_spacings.shape))
    temp_xlmatrix_to_xl_y = lambda xl_array: (xl_array * np.ones((ny,nxl))).transpose()
    
    ztop_xl_y = temp_xlmatrix_to_xl_y(ztop)
    z_new = temp_ycol_to_xl_y(std_y_positions)
    new_y_spacings = temp_ycol_to_xl_y(std_y_spacings)
    std_layer_bottoms_xl_y = temp_ycol_to_xl_y(std_layer_bottoms)
    z_new[std_layer_bottoms_xl_y >= ztop_xl_y] = np.nan
    new_y_spacings[np.isnan(z_new)] = np.nan
    
    top_layer_index = temp_ycol_to_xl_y(np.arange(ny))
    top_layer_index[np.isnan(z_new)] = np.nan
    top_layer_index = np.nanmin(top_layer_index,axis=1)
    tlis = (np.arange(nxl), top_layer_index.astype(int))
    new_y_spacings[tlis[0],tlis[1]] = ztop_xl_y[tlis[0],tlis[1]] - std_layer_bottoms_xl_y[tlis[0],tlis[1]]
    z_new[tlis[0],tlis[1]] = std_layer_bottoms_xl_y[tlis[0],tlis[1]] + new_y_spacings[tlis[0],tlis[1]] / 2.0
    
    dy[x_indices,:,l_indices] = new_y_spacings
    z_depth[x_indices,:,l_indices] = z_new
    
    sthuf_temp = sthu[x_indices,:,l_indices] + sthf[x_indices,:,l_indices]
    
    # Going straight for interpolate sthu
    
    soilvar_array = np.array([v_sat[x_indices,:,l_indices], bexp[x_indices,:,l_indices], sathh[x_indices,:,l_indices], hcap[x_indices,:,l_indices], hcon[x_indices,:,l_indices], xice[x_indices,:,l_indices], sthuf_temp, tsl[x_indices,:,l_indices], sthu[x_indices,:,l_indices], sthf[x_indices,:,l_indices]]) # Check all needed!
    
    ### --- Interpolate --- ###
    #not_isnan = ~np.isnan(z_new) # Hmm, are there any nans here? What does this array look like? Weirdly nans are important, where do they come from?!
    soilvar_array = interp_soilvar_multiarray(soilvar_array, z_eff, z_new, tlis)
    
    v_sat[x_indices,:,l_indices], bexp[x_indices,:,l_indices], sathh[x_indices,:,l_indices], hcap[x_indices,:,l_indices], hcon[x_indices,:,l_indices], xice[x_indices,:,l_indices], sthuf_temp, tsl_temp, sthu_temp, sthf_temp = soilvar_array
    
    # Correction - Limit xice
    xice[x_indices,:,l_indices][xice[x_indices,:,l_indices] < 0] = 0.0 # should check other vars at some point too 
    #-----------------------------
    # This bit is like well dodgy! I don't know, interpolating temperature, what was he thinking... (probably should do enthalpy instead)!
    # ice_time_check('g',then,check_xice_times)
    
    smclsat[x_indices,:,l_indices] = rho_water * dx[x_indices,:,l_indices] * (1 - xice[x_indices,:,l_indices])*dy[x_indices,:,l_indices] * v_sat[x_indices,:,l_indices]
    # Correction, limit sthuf - # note, maybe this should just be reset by driving as in surface temperatures. Check if other interpolating variables are not volumetric
    sthuf_temp[sthuf_temp>1.0] = 1.0
    smcl[x_indices,:,l_indices] = smclsat[x_indices,:,l_indices]*sthuf_temp
    
    tmax[x_indices,:,l_indices] = calc_tmax(v_sat[x_indices,:,l_indices], smcl[x_indices,:,l_indices], smclsat[x_indices,:,l_indices], bexp[x_indices,:,l_indices], sathh[x_indices,:,l_indices]) # only updating smcl monthly, note calc_tmax is array shape independent
    
    # Again, straight for interpolate sthu
    
    T_temp = -(sathh[x_indices,:,l_indices]/dpsidt) * sthu_temp ** (-bexp[x_indices,:,l_indices])
    # First attempt, see if there is any sthf about, if so, use this t_temp, otherwise use interpolated T
    # Maybe some thinking about limits to be done!
    # need equal to conditions
    
    temp_filter = sthf_temp >= 0.001 # frozen pore ice
    tsl[x_indices,:,l_indices][temp_filter] = T_temp[temp_filter]
    sthu[x_indices,:,l_indices][temp_filter] = sthu_temp[temp_filter]
    sthf[x_indices,:,l_indices][temp_filter] = sthuf_temp[temp_filter] - sthu[x_indices,:,l_indices][temp_filter]
    
    temp_filter = (sthf_temp < 0.001) * (xice[x_indices,:,l_indices] > small_value) # unfrozen pore water, but with xice
    tsl[x_indices,:,l_indices][temp_filter] = tmax[x_indices,:,l_indices][temp_filter] 
    sthf[x_indices,:,l_indices][temp_filter] = 0.0
    sthu[x_indices,:,l_indices][temp_filter] = sthuf_temp[temp_filter]
    
    
    temp_filter = (sthf_temp < 0.001) * (xice[x_indices,:,l_indices] <= small_value) # all unfrozen
    tsl[x_indices,:,l_indices][temp_filter] = tsl_temp[temp_filter]
    sthf[x_indices,:,l_indices][temp_filter] = 0.0
    sthu[x_indices,:,l_indices][temp_filter] = sthuf_temp[temp_filter]
    
    # End of interpolate sthu bit
    
    if no_heat_flux_edges:
        for x in [1, nx-2]:
            les = l_indices[x_indices == x]
            xe = {1:0, nx - 2: nx-1}[x]
            v_sat[xe,:,les] = v_sat[x,:,les]
            smclsat[xe,:,les] = smclsat[x,:,les]
            bexp[xe,:,les] = bexp[x,:,les]
            sathh[xe,:,les] = sathh[x,:,les]
            hcap[xe,:,les] = hcap[x,:,les]
            hcon[xe,:,les] = hcon[x,:,les]
            xice[xe,:,les] = xice[x,:,les]
            smcl[xe,:,les] = smcl[x,:,les]
            sthu[xe,:,les] = sthu[x,:,les] # smclu and smclf not calculated from these next timestep
            sthf[xe,:,les] = sthf[x,:,les]
            tsl[xe,:,les] = tsl[x,:,les]
            dy[xe,:,les] = dy[x,:,les]
            z_depth[xe,:,les] = z_depth[x,:,les]

    return xice, dy, z_depth, tsl,  tmax, smcl, smclsat, sthu, sthf, v_sat, bexp, sathh, hcap, hcon

def excess_ice_calc_original(xice, thawing_xice, dxice_tot, dy, z_depth,
                    nx, ny, nl, no_heat_flux_edges, dx,
                    std_y_positions, std_y_spacings, std_layer_tops,
                    tsl, tmax, smcl, smclu, smclf, smclsat, sthu, sthf, v_sat, bexp, sathh, hcap, hcon):
    """
    The one to use for the moment! Sorts out subsidence if xice has thawed.
    
    """
    check_consistency = True # For now

    for x, l in zip(np.indices((nx,nl))[0][thawing_xice.max(axis=1)],
                    np.indices((nx,nl))[1][thawing_xice.max(axis=1)]): # Could we do better than a for loop?
        if not((x == 0 or x == nx-1) and no_heat_flux_edges):
            # For speed, might not want to keep reinitialising this array
            sthuf_temp = sthu[x,:,l] + sthf[x,:,l]
            
            if check_consistency:
                init_array = np.array([v_sat[x,:,l], bexp[x,:,l], sathh[x,:,l], hcap[x,:,l], hcon[x,:,l], xice[x,:,l], sthuf_temp, tsl[x,:,l], dy[x,:,l], z_depth[x,:,l]])
            # Remember, need to define y_spacings - probably use dy?
            
            eff_y_spacings = dy[x,:,l] * (1 + dxice_tot[x,:,l])
            z_eff = tot_depth + thickness_to_positions(eff_y_spacings[::-1])[::-1]
            ztop = tot_depth + np.nansum(eff_y_spacings) # have we defined tot_depth?

            z_new = std_y_positions.copy()
            new_y_spacings = std_y_spacings.copy()
            z_new[std_layer_bottoms + small_value * ny >= ztop] = np.nan
            new_y_spacings[std_layer_bottoms + small_value * ny >= ztop] = np.nan
            z_new[np.isnan(z_depth[x,:,l])] = np.nan # Only one way at the moment
            new_y_spacings[np.isnan(z_depth[x,:,l])] = np.nan

            top_layer_index = np.min(np.arange(len(z_new))[~np.isnan(z_new)])
            new_y_spacings[top_layer_index] = ztop - std_layer_bottoms[top_layer_index]
            z_new[top_layer_index] = std_layer_bottoms[top_layer_index] + new_y_spacings[top_layer_index] / 2.0

            dy[x,:,l] = new_y_spacings
            z_depth[x,:,l] = z_new


            if not interpolate_sthu:
                soilvar_array = np.array([v_sat[x,:,l], bexp[x,:,l], sathh[x,:,l], hcap[x,:,l], hcon[x,:,l], xice[x,:,l], sthuf_temp, tsl[x,:,l]])
                soilvar_array_temp = np.full(np.shape(soilvar_array),np.nan)

                ### --- Interpolate --- ###
                if first_order_interp:
                    soilvar_array_temp[:,~np.isnan(z_new)] = interp_soilvar_1st_order(np.array([soilvar[~np.isnan(z_new)] for soilvar in soilvar_array]), z_eff[~np.isnan(z_new)], z_new[~np.isnan(z_new)])
                else:
                    soilvar_array_temp[:,~np.isnan(z_new)] = interp_soilvar_faster(np.array([soilvar[~np.isnan(z_new)] for soilvar in soilvar_array]), z_eff[~np.isnan(z_new)], z_new[~np.isnan(z_new)])
                v_sat[x,:,l], bexp[x,:,l], sathh[x,:,l], hcap[x,:,l], hcon[x,:,l], xice[x,:,l], sthuf_temp, tsl[x,:,l] = soilvar_array_temp

            else: # interpolate sthu
                soilvar_array = np.array([v_sat[x,:,l], bexp[x,:,l], sathh[x,:,l], hcap[x,:,l], hcon[x,:,l], xice[x,:,l], sthuf_temp, tsl[x,:,l], sthu[x,:,l], sthf[x,:,l]])
                soilvar_array_temp = np.full(np.shape(soilvar_array),np.nan)
    
                ### --- Interpolate --- ###
                if first_order_interp:
                    soilvar_array_temp[:,~np.isnan(z_new)] = interp_soilvar_1st_order(np.array([soilvar[~np.isnan(z_new)] for soilvar in soilvar_array]), z_eff[~np.isnan(z_new)], z_new[~np.isnan(z_new)])
                else:
                    soilvar_array_temp[:,~np.isnan(z_new)] = interp_soilvar_faster(np.array([soilvar[~np.isnan(z_new)] for soilvar in soilvar_array]), z_eff[~np.isnan(z_new)], z_new[~np.isnan(z_new)])
                v_sat[x,:,l], bexp[x,:,l], sathh[x,:,l], hcap[x,:,l], hcon[x,:,l], xice[x,:,l], sthuf_temp, tsl_temp, sthu_temp, sthf_temp = soilvar_array_temp
            
            # Correction - Limit xice
            xice[x,:,l][xice[x,:,l] < 0.0] = 0.0 # should check other vars at some point too 
            if correct_thaw:
                z_top = z_depth[x, top_layer_index, l] + dy[x, top_layer_index, l]/2.0
                z_error = z_top - np.nansum(dy[x,:,l]*xice[x,:,l])
                if np.sqrt(z_error**2) > small_value * ny:
                    if z_error < dy[x, top_layer_index, l] and dy[x, top_layer_index, l] - z_error < std_y_spacings[top_layer_index]:
                        dy[x, top_layer_index, l] -= z_error
                        z_depth[x, top_layer_index, l] = std_layer_bottoms[top_layer_index] + dy[x, top_layer_index, l] / 2.0
                        # For now, not removing or adding a layer, but ignoring it until this is possible?

            #-----------------------------
            # This bit is like well dodgy! I don't know, interpolating temperature, what was he thinking... (probably should do enthalpy instead)!
            
            smclsat[x,:,l] = rho_water * dx[x,:,l] * (1 - xice[x,:,l])*dy[x,:,l] * v_sat[x,:,l]
            # Correction, limit sthuf - # note, maybe this should just be reset by driving as in surface temperatures. Check if other interpolating variables are not volumetric
            sthuf_temp[sthuf_temp>1.0] = 1.0
            smcl[x,:,l] = smclsat[x,:,l]*sthuf_temp
            tmax[x,:,l] = calc_tmax(v_sat[x,:,l], smcl[x,:,l], smclsat[x,:,l], bexp[x,:,l], sathh[x,:,l]) # only updating smcl monthly
            
            if not interpolate_sthu:
                # Correction - limit soil temp while xice thawing
                to_limit = (xice[x,:,l] > 0) * (tsl[x,:,l] > tmax[x,:,l])
                tsl[x,:,l][to_limit] = tmax[x,:,l][to_limit] # see what happens
                
                sel_from_col = tsl[x,:,l] >= tmax[x,:,l]
                smclu[x,:,l][sel_from_col] = smcl[x,:,l][sel_from_col] # <!C22D!> # Unfrozen moisture content of each cell (kg/m).
                smclf[x,:,l][sel_from_col] = 0.0 # <!C22D!> # Frozen moisture content of each cell (kg/m).

                sel_from_col = tsl[x,:,l] < tmax[x,:,l]
                smclu[x,:,l][sel_from_col] = smclsat[x,:,l][sel_from_col] * (-dpsidt * tsl[x,:,l][sel_from_col] / sathh[x,:,l][sel_from_col])**(-1.0 / bexp[x,:,l][sel_from_col])
                smclf[x,:,l][sel_from_col] = smcl[x,:,l][sel_from_col] - smclu[x,:,l][sel_from_col] # Hmm, what if smcl < smclu?!

                sthu  = smclu / smclsat
                sthf  = smclf / smclsat
                
            else:
                T_temp = -(sathh[x,:,l]/dpsidt) * sthu_temp ** (-bexp[x,:,l])
                # First attempt, see if there is any sthf about, if so, use this t_temp, otherwise use interpolated T
                # Maybe some thinking about limits to be done!

                temp_filter = sthf_temp >= 0.001 # frozen pore ice
                tsl[x,:,l][temp_filter] = T_temp[temp_filter]
                sthu[x,:,l][temp_filter] = sthu_temp[temp_filter]
                sthf[x,:,l][temp_filter] = sthuf_temp[temp_filter] - sthu[x,:,l][temp_filter]
                
                temp_filter = (sthf_temp < 0.001) * (xice[x,:,l] > small_value) # unfrozen pore water, but with xice
                tsl[x,:,l][temp_filter] = tmax[x,:,l][temp_filter] 
                sthf[x,:,l][temp_filter] = 0.0
                sthu[x,:,l][temp_filter] = sthuf_temp[temp_filter]

                
                temp_filter = (sthf_temp < 0.001) * (xice[x,:,l] <= small_value) # all unfrozen
                tsl[x,:,l][temp_filter] = tsl_temp[temp_filter]
                sthf[x,:,l][temp_filter] = 0.0
                sthu[x,:,l][temp_filter] = sthuf_temp[temp_filter]

                tsl[x,:,l][np.isnan(tsl_temp)] = np.nan
                sthu[x,:,l][np.isnan(sthu_temp)] = np.nan # could be quicker
                sthf[x,:,l][np.isnan(sthf_temp)] = np.nan # The consistency checks might not work now.
            
            if (x == 1 or x == nx - 2) and no_heat_flux_edges:
                xe = {1:0, nx - 2: nx-1}[x]
                v_sat[xe,:,l] = v_sat[x,:,l]
                smclsat[xe,:,l] = smclsat[x,:,l]
                bexp[xe,:,l] = bexp[x,:,l]
                sathh[xe,:,l] = sathh[x,:,l]
                hcap[xe,:,l] = hcap[x,:,l]
                hcon[xe,:,l] = hcon[x,:,l]
                xice[xe,:,l] = xice[x,:,l]
                smcl[xe,:,l] = smcl[x,:,l]
                sthu[xe,:,l] = sthu[x,:,l] # smclu and smclf not calculated from these next timestep
                sthf[xe,:,l] = sthf[x,:,l]
                tsl[xe,:,l] = tsl[x,:,l]
                dy[xe,:,l] = dy[x,:,l]
                z_depth[xe,:,l] = z_depth[x,:,l]

            if check_consistency:
                final_array = np.array([v_sat[x,:,l], bexp[x,:,l], sathh[x,:,l], hcap[x,:,l], hcon[x,:,l], xice[x,:,l], sthuf_temp, dy[x,:,l], z_depth[x,:,l]])
                checks = np.array([(np.isnan(tsl[x,:,l]) != np.isnan(soilvar)).any() for soilvar in final_array])
                if checks.any():
                    print(f'Consistency checks failed for x{x} and l{l}!')
                    print(f"Nanpattern for tsl != nanpattern for {np.array(['v_sat', 'bexp', 'sathh', 'hcap', 'hcon', 'xice', 'sthuf_temp', 'dy', 'z_depth'])[checks]}")
                    print('Initial array:')
                    for varname, var in zip(init_array,['v_sat', 'bexp', 'sathh', 'hcap', 'hcon', 'xice', 'sthuf_temp', 'tsl','dy', 'z_depth']):
                        print(var, varname)
                    print('Final array:')
                    for varname, var in zip(final_array,['v_sat', 'bexp', 'sathh', 'hcap', 'hcon', 'xice', 'sthuf_temp', 'dy', 'z_depth']):
                        print(var, varname)
                    print('tsl:')
                    print(tsl[x,:,l])
                    raise Exception('This ends here (new).')

    return xice, dy, z_depth, tsl, tmax, smcl, smclsat, sthu, sthf, v_sat, bexp, sathh, hcap, hcon
