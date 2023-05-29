"""
Global 2D

"""
code_test = False
global_test = False

print('global_2D_current last updated 28/01/2023 22:34')

# Useful for debugging:
# if current_time == datetime.datetime(1901,10,1,0,0):
#     print('sleep...')
#     time.sleep(30)
        
#----------------------------------
# Parallel input
#----------------------------------
param_option = 'default'
subchunk = False
if not code_test:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--point', nargs=2,type=float,
                        help='Runs the lat (degrees N) and lon (degrees E) specified.')
    parser.add_argument('--chunk', nargs=2,type=int, 
                        help='Runs the chunk specified by the index of task and number of tasks')
    parser.add_argument('--point_option', nargs=1,type=int, 
                        help='Specifies the set of parameter options.')
    parser.add_argument('--subchunk', nargs=4,type=int, 
                        help='Runs a subchunk of the chunk specified by the old index of task, old number of tasks, thesubchunk index and number of subchunks')
    args = parser.parse_args()
    if args.point_option != None:
        point_run = True
        point_lat, point_lon, param_option = args.point[0], args.point[1], args.point_option[0]
        print(f'Point Option run. Point for run: {point_lat}N, {point_lon}E, param_option {param_option}')
        task_n, n_tasks, subtask_n, n_subtasks = None, None, None, None
    elif args.point != None:
        point_run = True
        point_lat, point_lon = args.point[0], args.point[1]
        print(f'Point run. Point for run: {point_lat}N, {point_lon}E')
        task_n, n_tasks, subtask_n, n_subtasks = None, None, None, None
    elif args.chunk != None:
        point_run = False
        task_n, n_tasks, subtask_n, n_subtasks = args.chunk[0], args.chunk[1], None, None
        print(f'Chunk run. Not a point. Task number {task_n} of {n_tasks} tasks.')
        point_lat, point_lon = None, None
    elif args.subchunk != None:
        point_run = False
        subchunk = True
        task_n, n_tasks, subtask_n, n_subtasks = args.subchunk[0], args.subchunk[1], args.subchunk[2], args.subchunk[3]
        print(f'Subchunk run. Task number {task_n} of {n_tasks} tasks, subtask {subtask_n} of {n_subtasks} subtasks.')
        point_lat, point_lon = None, None
else:
    if global_test:
        point_run = False
        subchunk = True
        task_n, n_tasks, subtask_n, n_subtasks = 192, 200, 2, 10 # for prefixes
        print(f'Global test. Cropping latlon.')
        point_lat, point_lon = None, None
    else:
        point_run = True
        point_lat, point_lon = 69.25, 25.25
        print(f'Point test. Point for run: {point_lat}N, {point_lon}E')
        task_n, n_tasks, subtask_n, n_subtasks = None, None, None, None

if not code_test:
    import pickle # Terrible, but do a classier approach some other time.
    from parm_options_mod import parm_options_dict
    parm_options_sel = parm_options_dict[param_option]
    if point_run: # ahhh!
        parm_options_sel['prefix_suffix'] = f'_{point_lat}N_{point_lon}E_{param_option}'
    else: # though this probably doesn't make sense anyway
        if subchunk:
            parm_options_sel['prefix_suffix'] = f'_{task_n:03}_{subtask_n:02}_{param_option}'
        else:
            parm_options_sel['prefix_suffix'] = f'_{task_n:03}_{param_option}'
    parm_pickle_file = open('sel_parm_options.pk', 'ab')
    pickle.dump(parm_options_sel, parm_pickle_file)
    parm_pickle_file.close()

#----------------------------------
# Import definitions and functions
#----------------------------------
if not code_test:
    from global_lat_thaw_parm_mod import *
import pprint
print(f"In global_2D_current, sp =")
pprint.pp(vars(sp))
if point_run:
    file_prefix = f'{prefix_prefix}_{point_lat}N_{point_lon}E_{sp.param_option}'
else:
    if subchunk:
        file_prefix = f'{prefix_prefix}_{task_n:03}_{subtask_n:02}_{sp.param_option}'
    else:
        file_prefix = f'{prefix_prefix}_{task_n:03}_{sp.param_option}'

print(f'File prefix = {file_prefix}')

#----------------------------------
# Initialise
#----------------------------------
print('Initialise_input')
nl, dx, dy, soil_depths, high, low, high_daily, low_daily, v_sat, bexp, sathh, hcap, hcon, smclsat, tsl, smcl, tmax, smclu, smclf, sthu, sthf, n_unsat, xice, z_depth, is_surface, T_surf_gradient_grid, T_surf_high_grid = initialise_input(task_n, n_tasks, subtask_n, n_subtasks, file_prefix, point_run, point_lat, point_lon, subchunk)

#start_datetime = datetime.datetime(2000,1,1,0,0)
#end_datetime = datetime.datetime(2010,1,1,0,0)
timestep_timedelta = datetime.timedelta(seconds=timestep)
def process_datetime(date):
    date = np.datetime64(date,'M')
    date = pd.Timestamp(date - np.timedelta64(1, 'M')).to_pydatetime()
    return date
start_datetime = process_datetime(high.time[0].values)
end_datetime = process_datetime(high.time[-1].values)
timestep_datetimes = np.arange(start_datetime, end_datetime, timestep_timedelta).astype(datetime.datetime)

print('Initialise output')
outputs, total_timedeltas  = initialise_outputs(output_profiles, end_datetime, 
                                start_datetime, soil_depths, nx, ny, nl, thickness, high)
print('Initisation complete.')

#-----------------------------------------------------------------------------------
# Main run
#-----------------------------------------------------------------------------------
print('Starting main run')
                                                           
last_save = 0
n_days = 0
n_months = 0
n_years = 0
day_start = time.time()
main_start = time.time()

daily_lat_heat = np.zeros((ny,nl)) # Now layered
monthly_lat_heat = np.zeros((ny,nl))
yearly_lat_heat = np.zeros((ny,nl))
daily_remained_frozen = np.zeros((nx,ny,nl)) # Records if remained frozen since last refreshed.
monthly_remained_frozen = np.zeros((nx,ny,nl))
yearly_remained_frozen = np.zeros((nx,ny,nl))
daily_tsl = np.zeros((nx,ny,nl))
monthly_tsl = np.zeros((nx,ny,nl))
yearly_tsl = np.zeros((nx,ny,nl))
x_offset = np.zeros(nl)
x_of_y_threshold_lower = np.zeros(nl)
x_of_y_threshold_upper = np.zeros(nl)
prev_i_d = 0
prev_i_m = 0
prev_i_y = 0

temp_grid = np.indices((ny,nl))
yl_y = temp_grid[0]
yl_l = temp_grid[1]
temp_grid = np.indices((nx,nl))
xl_x = temp_grid[0]
xl_l = temp_grid[1]
previous_time = timestep_datetimes[0]
# with tqdm, smoothing = 0 for an average,1 for instantaneous, default is 0.3
if time_to_go == 'tqdm':
    iterator = enumerate(tqdm(timestep_datetimes, smoothing = 0))
else:
    iterator = enumerate(timestep_datetimes)

for i, current_time in iterator: # or enumerate(tqdm(timestep...
    #-------------------------------
    # Maybe update driving
    #-------------------------------
    if previous_time != current_time:
        tsl, smcl, tmax, T_surf_gradient_grid, T_surf_high_grid = update_driving(previous_time, current_time, tsl, smcl, high, low, 
                                         smclsat, tmax, high_daily, low_daily, n_unsat, v_sat, bexp, sathh, z_depth, is_surface, T_surf_gradient_grid,T_surf_high_grid, nl)
    if no_heat_flux_edges: # could do something neater...
        tsl[0,:,:] = tsl[1,:,:]
        tsl[-1,:,:] = tsl[-2,:,:]
        tsl[:,-1,:] = tsl[:,-2,:]
    if sp.option == 'twotilesat':
        tsl[:,-1,:] = tsl[:,-2,:]
    smclu = smclsat * sthu # Frozen moisture content of each soil layer (kg/m2).
    smclf = smcl - smclu # Have changed this round relative to JULES!
    
    #-------------------------------------------------------------------------
    # Calculate conductivities and heat flux to active points
    #-------------------------------------------------------------------------
    if sp.option == 'twotilesat': # no x
        #swapped indicies to be like jules - should it be like this?
        smcfu_y=(dy[:,1:,:] * sthu[:,:-1,:] + dy[:,:-1,:] * sthu[:,1:,:]) / (dy[:,:-1,:] + dy[:,1:,:])
        smcff_y=(dy[:,1:,:] * sthf[:,:-1,:] + dy[:,:-1,:] * sthf[:,1:,:]) / (dy[:,:-1,:] + dy[:,1:,:])
        # like jules - for bit comparability, but why?
        #v_satk_y=(dy[:,:-1,:] * v_sat[:,:-1,:] + dy[:,1:,:] * v_sat[:,1:,:]) / (dy[:,:-1,:] + dy[:,1:,:])
        # v_sat currently used in JULES for layer boundary, but not sure it should be?
        v_satk_y = v_sat[:,:-1,:]
        hcon_y = hcon[:,:-1,:]
        # as the hcon is already for between layers, 
        # but why is it the same length as the number of soil layers? 
        # -> last one for link to deep layer. 
        k_y_grid = heat_con(hcon_y, smcfu_y, smcff_y, v_satk_y)

        if sp.excess_ice:
            # Assuming horizontal ice lenses, series in y so resistances add
            avg_xice = (dy[:,1:,:] * xice[:,:-1,:] + dy[:,:-1,:] * xice[:,1:,:]) / (dy[:,:-1,:] + dy[:,1:,:])
            k_y_grid = 1/((1 - avg_xice) * (1/ k_y_grid) + avg_xice * (1 / hcice))

        dhsl = np.zeros(shape=(nx,ny,nl))
        dhsl[:,1:-1,:] = ( k_y_grid[:,1:,:]*dx[:,1:-1,:]*(tsl[:,2:,:] - tsl[:,1:-1,:]) / (0.5*(dy[:,2:,:] + dy[:,1:-1,:]))
                           -k_y_grid[:,:-1,:]*dx[:,1:-1,:]*(tsl[:,1:-1,:] - tsl[:,:-2,:]) / (0.5*(dy[:,1:-1,:] + dy[:,:-2,:])) )
    else:
        #swapped indicies to be like jules - should it be like this?
        smcfu_x=(dx[1:,:,:] * sthu[:-1,:,:] + dx[:-1,:,:] * sthu[1:,:,:]) / (dx[:-1,:,:] + dx[1:,:,:])
        # Fractional saturation (unfrozen water) at cell boundaries in x direction.
        smcff_x=(dx[1:,:,:] * sthf[:-1,:,:] + dx[:-1,:,:] * sthf[1:,:,:]) / (dx[:-1,:,:] + dx[1:,:,:])
        #v_satk_x=(dx[:-1,:,:] * v_sat[:-1,:,:] + dx[1:,:,:] * v_sat[1:,:,:]) / (dx[:-1,:,:] + dx[1:,:,:])
        #hcon_x=(dx[:-1,:,:] * hcon[:-1,:,:] + dx[1:,:,:] * hcon[1:,:,:]) / (dx[:-1,:,:] + dx[1:,:,:])
        v_satk_x = v_sat[:-1,:,:] # see below comments, probably fine for x since uniform soil properties.
        hcon_x = hcon[:-1,:,:] # however, if the hcon is for layer bounds, maybe this ought be an average of the upper and lower bounds?
        k_x_grid = heat_con(hcon_x, smcfu_x, smcff_x, v_satk_x)

        if sp.excess_ice:
            # Assuming horizontal ice lenses, so conductivities in parallel in x so add
            # Note, just taking average excess ice though.
            avg_xice = (dx[1:,:,:] * xice[:-1,:,:] + dx[:-1,:,:] * xice[1:,:,:]) / (dx[:-1,:,:] + dx[1:,:,:])
            k_x_grid = hcice * avg_xice + k_x_grid * (1 - avg_xice)
        
        #swapped indicies to be like jules - should it be like this?
        smcfu_y=(dy[:,1:,:] * sthu[:,:-1,:] + dy[:,:-1,:] * sthu[:,1:,:]) / (dy[:,:-1,:] + dy[:,1:,:])
        smcff_y=(dy[:,1:,:] * sthf[:,:-1,:] + dy[:,:-1,:] * sthf[:,1:,:]) / (dy[:,:-1,:] + dy[:,1:,:])
#         v_satk_y=(dy[:,:-1,:] * v_sat[:,:-1,:] + dy[:,1:,:] * v_sat[:,1:,:]) / (dy[:,:-1,:] + dy[:,1:,:])
        # v_sat currently used in JULES for layer boundary, but not sure it should be?
        v_satk_y = v_sat[:,:-1,:]
        hcon_y = hcon[:,:-1,:]
        # as the hcon is already for between layers, 
        # but why is it the same length as the number of soil layers? 
        # -> last one for link to deep layer. 
        k_y_grid = heat_con(hcon_y, smcfu_y, smcff_y, v_satk_y) # ALL the conductivities in the x direction

        if sp.excess_ice:
            # Assuming horizontal ice lenses, series in y so resistances add
            avg_xice = (dy[:,1:,:] * xice[:,:-1,:] + dy[:,:-1,:] * xice[:,1:,:]) / (dy[:,:-1,:] + dy[:,1:,:])
            k_y_grid = 1/((1 - avg_xice) * (1/ k_y_grid) + avg_xice * (1 / hcice))

        dhsl = np.zeros(shape=(nx,ny,nl))
        dhsl[1:-1,1:-1,:] = ( k_x_grid[1:,1:-1,:]*dy[1:-1,1:-1,:]*(tsl[2:,1:-1,:] - tsl[1:-1,1:-1,:]) / (0.5*(dx[2:,1:-1,:] + dx[1:-1,1:-1,:]))
                           -k_x_grid[:-1,1:-1,:]*dy[1:-1,1:-1,:]*(tsl[1:-1,1:-1,:] - tsl[:-2,1:-1,:]) / (0.5*(dx[1:-1,1:-1,:] + dx[:-2,1:-1,:]))
                           +k_y_grid[1:-1,1:,:]*dx[1:-1,1:-1,:]*(tsl[1:-1,2:,:] - tsl[1:-1,1:-1,:]) / (0.5*(dy[1:-1,2:,:] + dy[1:-1,1:-1,:]))
                           -k_y_grid[1:-1,:-1,:]*dx[1:-1,1:-1,:]*(tsl[1:-1,1:-1,:] - tsl[1:-1,:-2,:]) / (0.5*(dy[1:-1,1:-1,:] + dy[1:-1,:-2,:])) )
        #if option == 'offset': dhsl[np.arange(nx),surf_y_index] = 0.0
        dhsl[np.isnan(dhsl)] = 0.0 # might work
        dhsl[is_surface] = 0.0 # just to make sure?
        # The heat available to update the cell temperature (J/m/timestep)
        # Here we choose not to have heat flow between outer edge cells.
        
        #Lat thaw calc
        here = np.indices((nx,ny,nl))[0]
        here[tsl>=0] = 0 # did have nx and used a min
        here[np.isnan(tsl)]=0
        i_x = here.max(axis=0).astype(int) # horizontal position of thaw front
        #i_x[i_x==0] = np.nan # mire frozen, could use a distinguishing value?
        i_x[i_x==nx-1] = 0 # all unfrozen, could use a distinguishing value?
        dhslx = np.zeros((nx,ny,nl)) # could make more efficient
        dhslx[1:-1,1:-1,:] = ( k_x_grid[1:,1:-1,:]*dy[1:-1,1:-1,:]*(tsl[2:,1:-1,:] - tsl[1:-1,1:-1,:]) / (0.5*(dx[2:,1:-1,:] + dx[1:-1,1:-1,:]))
                           -k_x_grid[:-1,1:-1,:]*dy[1:-1,1:-1,:]*(tsl[1:-1,1:-1,:] - tsl[:-2,1:-1,:]) / (0.5*(dx[1:-1,1:-1,:] + dx[:-2,1:-1,:])) )
        lat_heat_layers = dhsl[i_x[:,:],yl_y,yl_l]
        daily_lat_heat = daily_lat_heat + timestep*lat_heat_layers

    # Skipped heat flux bit for now

    points, thawing_xice, itmax = np.full((nx,ny,nl), True), np.full((nx,ny,nl), False), 3
    dxice_tot = np.zeros(shape=(nx,ny,nl)) # for iteration. could make more efficient
    original_xice = xice.copy() # For use in heat cpacity, so that the amount of soil doesn't grow
    
    while points.any() and itmax != 0: # Iteration could be more efficient as lots of double calculating.
        #--------------------------------
        # Phase changes and heat capacity
        #--------------------------------
        unfrozen = tsl > tmax
        frozen = tsl < tmax

        # dthu = rate of change of volumetric unfrozen soil moisture concentration with temperature (m3 liquid H2O/m3 soil/K)
        dthu = np.zeros(shape=(nx,ny,nl))
        calc_dthu = frozen + (np.sqrt((tsl-tmax)**2) < small_value)*(dhsl < 0.0) # Frozen or on the edge and freezing # Rounding 1
        dthu[calc_dthu] = ( dpsidt * v_sat[calc_dthu] / (bexp[calc_dthu] * sathh[calc_dthu]) 
                      * (-dpsidt * tsl[calc_dthu] / sathh[calc_dthu])**(-1.0 / bexp[calc_dthu] - 1.0) ) # Frozen bannanas

        if sp.excess_ice:
            xice_thaw = (tsl >= (tmax - small_value)) * (dhsl >= 0.0) * (xice > small_value) # Rounding 2 - do these small_value(s) work?
            if no_heat_flux_edges:
                xice_thaw[[0,-1],:,:] = False # Speed, avoid calculating. Also, dhsl should be 0.
                xice_thaw[:,-1,:] = False
            dxice_temp = - np.fmin( (dhsl[xice_thaw] * timestep / 
                        (lf * rho_water * v_sat[xice_thaw] * dx[xice_thaw] * dy[xice_thaw]) ), 
                                xice[xice_thaw])
            xice[xice_thaw] += dxice_temp
            dhsl[xice_thaw] -= ( dxice_temp * lf * rho_water * v_sat[xice_thaw] * 
                                dx[xice_thaw] * dy[xice_thaw] / timestep )
            thawing_xice[xice_thaw] = True
            dxice_tot[xice_thaw] += dxice_temp

        # hcapt = total volumetric heat capacity (soil+water) of the layer (J/m3/K).
        hcapt = ((1 - original_xice)*( hcap + (hcapw * sthu + hcapi * sthf) * rho_water * v_sat # Hmm, I've changed this, but I wonder if it was important to take into account that ice is less dense?
                  + rho_water * dthu * ((hcapw - hcapi) * tsl + lf) )
                  + hcapi * rho_water * xice )

        #--------------------------------
        # Add heat
        #--------------------------------
        dtsl0 = timestep*dhsl/(hcapt*dx*dy)
        dtsl = dtsl0.copy()

        dtslim = tmax - tsl
        
        temp_points_unfrozen = np.full((nx,ny,nl), False)
        temp_points_frozen = np.full((nx,ny,nl), False)
        temp_points_unfrozen[unfrozen] = dtsl[unfrozen] < dtslim[unfrozen] # Unfrozen, but about to freeze
        temp_points_frozen[frozen] = dtsl[frozen] > dtslim[frozen] # Frozen, but about to thaw
        
        dtsl[temp_points_unfrozen] = dtslim[temp_points_unfrozen] # could do double index?
        dtsl[temp_points_frozen]   = dtslim[temp_points_frozen] # Limiting the temperature change for this iteration
        
        dhsl = (dtsl0 - dtsl)*(hcapt*dx*dy)/timestep # Heat carried forward to the next iteration
        
        # could cut below for speed?
        if (dtsl > (100.0 - tsl)).any(): raise Exception('Boiling!')
        if (dtsl < (- zerodegc - tsl)).any(): raise Exception('Absolute 0!')

        tsl[points] = tsl[points] + dtsl[points] # not sure if necessary to filter by points really.
        points = np.logical_and(points, np.logical_or(temp_points_unfrozen, temp_points_frozen))
        
        # Skipped what to do if small temp increment and frozen water exists.

        #-----------------------------------------------------------------------
        # Diagnose unfrozen and frozen water contents
        #-----------------------------------------------------------------------

        smclu[tsl >= tmax] = smcl[tsl >= tmax] # <!C22D!> # Unfrozen moisture content of each cell (kg/m).
        smclf[tsl >= tmax] = 0.0 # <!C22D!> # Frozen moisture content of each cell (kg/m).

        smclu[tsl<tmax] = smclsat[tsl<tmax] * (-dpsidt * tsl[tsl<tmax] / sathh[tsl<tmax])**(-1.0 / bexp[tsl<tmax])
        smclf[tsl<tmax] = smcl[tsl<tmax] - smclu[tsl<tmax] # Hmm, what if smcl < smclu?!

        # Skip check error for now.

        sthu  = smclu / smclsat
        sthf  = smclf / smclsat

        itmax -= 1
        
    # check excess ice thaw - maybe should be a separate function?
    if sp.excess_ice and thawing_xice.any():
        if excess_ice_code == 'multiarray':
            xice, dy, z_depth, tsl, tmax, smcl, smclsat, sthu, sthf, v_sat, bexp, sathh, hcap, hcon = excess_ice_calc_multiarray(xice, thawing_xice, dxice_tot, dy, z_depth,
                    nx, ny, nl, no_heat_flux_edges, dx,
                    std_y_positions, std_y_spacings, std_layer_tops,
                    tsl, tmax, smcl, smclu, smclf, smclsat, sthu, sthf, v_sat, bexp, sathh, hcap, hcon)

        elif excess_ice_code == 'original':
            xice, dy, z_depth, tsl, tmax, smcl, smclsat, sthu, sthf, v_sat, bexp, sathh, hcap, hcon = excess_ice_calc_original(xice, thawing_xice, dxice_tot, dy, z_depth,
                    nx, ny, nl, no_heat_flux_edges, dx,
                    std_y_positions, std_y_spacings, std_layer_tops,
                    tsl, tmax, smcl, smclu, smclf, smclsat, sthu, sthf, v_sat, bexp, sathh, hcap, hcon)
        else:
            raise Exception(f'{excess_ice_code} not an option for excess_ice_code')

        if thawing_xice.any():
            is_surface = find_surface_cells(z_depth) # Could make faster by only checking changes?

            if sp.option == 'offset' and sp.excess_ice: # Change so not doing every time in future? Also, maybe should change smcl?
                tsl[is_surface] = T_surf_gradient_grid[is_surface]*(sp.high_bit_elevation - std_y_spacings[0]/2 - z_depth[is_surface]) / sp.high_bit_elevation + T_surf_high_grid[is_surface] # Assuming surface layer depths the same in driving and here
                smcl[:,n_unsat:,:] = smclsat[:,n_unsat:,:] # will this be consistent? e.g. we were interpolating sthu, but we calculate sthu based on sthf in the next cycle, might this change sthu and so energy and temperature? - have changed round so keeps track of sthu at the beginning. Could just do changes?

    #-------------------------------------------------------------
    # Now saving stuff after it's updated, and incrementing dates.
    #-------------------------------------------------------------
    # Note: hcapt output may be wrong due to iteration.
    
    daily_remained_frozen += (tsl <= 0.0) # does this need to be the thaw max thing?
    daily_tsl += tsl
    
    if previous_time != current_time: # Will break for a timestep >= daily
        add_day = current_time.day != previous_time.day
        add_month = current_time.month != previous_time.month
        if add_day: # Maybe record stuff
            times = ['daily']
            monthly_tsl += daily_tsl
            yearly_tsl += daily_tsl
            monthly_lat_heat += daily_lat_heat
            yearly_lat_heat += daily_lat_heat
            monthly_remained_frozen += daily_remained_frozen
            yearly_remained_frozen += daily_remained_frozen
            
            daily_tsl_avg = daily_tsl/(i-prev_i_d)
            daily_remained_frozen_avg = daily_remained_frozen/(i-prev_i_d)
            monthly_remained_frozen_avg, yearly_remained_frozen_avg = None, None # Overwritten when needed.
            monthly_tsl_avg, yearly_tsl_avg = None, None
            if add_month:
                times.append('monthly')
                monthly_tsl_avg = monthly_tsl/(i-prev_i_m)
                monthly_remained_frozen_avg = monthly_remained_frozen/(i-prev_i_m)
            if current_time.year != previous_time.year:
                times.append('yearly')
                yearly_tsl_avg = yearly_tsl/(i-prev_i_y)
                yearly_remained_frozen_avg = yearly_remained_frozen/(i-prev_i_y)

            for what_time in times:
                time_n = {'daily':n_days,'monthly':n_months, 'yearly':n_years}[what_time]
                time_lat_heat = {'daily':daily_lat_heat,'monthly':monthly_lat_heat, 'yearly':yearly_lat_heat}[what_time]
                time_remained_frozen = {'daily':daily_remained_frozen_avg,'monthly':monthly_remained_frozen_avg, 'yearly':yearly_remained_frozen_avg}[what_time]
                time_tsl_avg = {'daily':daily_tsl_avg,'monthly':monthly_tsl_avg, 'yearly':yearly_tsl_avg}[what_time]
                for output in outputs:
                    if output['freq'] == what_time:
                        if output['vars'] =='ends':
                            xsclice = [0,-1]
                        else:
                            xsclice = slice(None, None)
                        z_depth_surface_calculated = False
                        lower_threshold_calculated = False
                        upper_threshold_calculated = False
                        for var in output['data'].data_vars:
                            # Special cases
                            if var == 'lat_thaw':
                                output['data'][var][time_n,:] = time_lat_heat.sum(axis=0) # hopefully the y axis, leaving shape(nl)
                            elif var == 'lat_thaw_layers':
                                output['data'][var][time_n,:,:] = time_lat_heat
                            elif var == 'remained_frozen':
                                output['data'][var][time_n,:,:] = time_remained_frozen
                            elif var == 'sthuf':
                                output['data'][var][time_n,:,:,:] = sthf[xsclice,:,:] + sthu[xsclice,:,:]
                            elif var == 'z_depth_surface':
                                here = np.indices((nx,ny,nl))[1]
                                here[np.isnan(tsl)] = 3000
                                i_y = here.min(axis=1).astype(int)
                                output['data'][var][time_n,:,:] = (z_depth[xl_x,i_y,xl_l] + dy[xl_x,i_y,xl_l]/2)[xsclice,:]
                            elif var == 'depth_unfrozen':
                                for k in [0, -1]: # probably could be faster
                                    here = np.indices((ny,nl))[0]
                                    here[tsl[k,:,:]>=0.0] = ny-1 # .compute() may be necessary for dask array boolean slicing
                                    i_depth = here.min(axis=0).astype(int)
                                    output['data'][var][time_n,k,:] = np.array([soil_depths[l] for l in i_depth])
                            # Everything else
                            elif var == 't_soil':
                                if avg_t_soil:
                                    output['data'][var][time_n,:,:,:] = time_tsl_avg[xsclice,:,:]
                                else:
                                    output['data'][var][time_n,:,:,:] = tsl[xsclice,:,:] # could rename tsl to t_soil so globals works?
                            elif var == 'x_of_y_threshold_lower' or var == 'x_of_y_threshold_upper' or var == 'z_depth_surface' or var == 'slope_gradient':
                                if not z_depth_surface_calculated:
                                    here = np.indices((nx,ny,nl))[1]
                                    here[np.isnan(tsl)] = 3000
                                    i_y = here.min(axis=1).astype(int)
                                    output['data']['z_depth_surface'][time_n,:,:] = (z_depth[xl_x,i_y,xl_l] + dy[xl_x,i_y,xl_l]/2)[xsclice,:]
                                    z_depth_surface_calculated = True
                                if var == 'x_of_y_threshold_lower' or (var == 'slope_gradient' and not lower_threshold_calculated):
                                    x_of_y_threshold_lower = np.array([interpolate.interp1d(output['data']['z_depth_surface'][time_n,::-1,n_l], std_x_positions[::-1], bounds_error=False, assume_sorted = True)(sp.y_threshold_lower) for n_l in range(nl)])
                                    output['data']['x_of_y_threshold_lower'][time_n,:] = x_of_y_threshold_lower
                                    lower_threshold_calculated = True
                                if var == 'x_of_y_threshold_upper' or (var == 'slope_gradient' and not upper_threshold_calculated):
                                    x_of_y_threshold_upper = np.array([interpolate.interp1d(output['data']['z_depth_surface'][time_n,::-1,n_l], std_x_positions[::-1], bounds_error=False, assume_sorted = True)(sp.y_threshold_upper) for n_l in range(nl)])
                                    output['data']['x_of_y_threshold_upper'][time_n,:] = x_of_y_threshold_upper
                                    upper_threshold_calculated = True
                                if var == 'slope_gradient': # everything should be calculated by now
                                    output['data'][var][time_n,:] = (sp.y_threshold_upper - sp.y_threshold_lower)/(x_of_y_threshold_upper - x_of_y_threshold_lower)
                            elif var == 'x_offset':
                                output['data'][var][time_n,:] = x_offset[:]
                            elif var == 'xice_total':
                                output['data'][var][time_n,:,:] = np.nansum(xice[xsclice,:,:]*dy[xsclice,:,:], axis=1)
                            else:
                                output['data'][var][time_n,:,:,:] = globals()[var][xsclice,:,:] # Slightly cheeky use of globals here.

            if 'daily' in times: 
                daily_lat_heat[:,:] = 0.0 # reset cumulative
                daily_remained_frozen = np.zeros((nx,ny,nl))
                daily_tsl[:,:] = 0.0
                prev_i_d = i
            if 'monthly' in times: 
                monthly_lat_heat[:,:] = 0.0
                monthly_remained_frozen = np.zeros((nx,ny,nl))
                monthly_tsl[:,:] = 0.0
                prev_i_m = i
            if 'yearly' in times: 
                yearly_lat_heat[:,:] = 0.0
                yearly_remained_frozen = np.zeros((nx,ny,nl))
                yearly_tsl[:,:] = 0.0
                prev_i_y = i
            if add_month: # Increment months, and maybe years, shift domain if needed.
                if sp.shifting_domain:
                    temp_x_spacing = dx[-1,-1,0] # assumes equal grid spacing
                    n_shift = int(2*(-sp.domain_x_threshold) / temp_x_spacing) # deciding to shift in chunks as will look better for plotting in 2D, but could do a single cell at a time?
                    for i_l in range(nl): # vectorise?
                        if x_of_y_threshold_lower[i_l] < sp.domain_x_threshold:
                            # This method is vulnerable to the code changing - another reason a soil grid object would be a good idea!
                            dx[n_shift:,:,i_l], dy[n_shift:,:,i_l], v_sat[n_shift:,:,i_l], bexp[n_shift:,:,i_l], sathh[n_shift:,:,i_l], hcap[n_shift:,:,i_l], hcon[n_shift:,:,i_l], smclsat[n_shift:,:,i_l], tsl[n_shift:,:,i_l], smcl[n_shift:,:,i_l], tmax[n_shift:,:,i_l], smclu[n_shift:,:,i_l], smclf[n_shift:,:,i_l], sthu[n_shift:,:,i_l], sthf[n_shift:,:,i_l], xice[n_shift:,:,i_l], z_depth[n_shift:,:,i_l], is_surface[n_shift:,:,i_l], daily_remained_frozen[n_shift:,:,i_l], monthly_remained_frozen[n_shift:,:,i_l], yearly_remained_frozen[n_shift:,:,i_l], daily_tsl[n_shift:,:,i_l], monthly_tsl[n_shift:,:,i_l], yearly_tsl[n_shift:,:,i_l] = dx[:-n_shift,:,i_l], dy[:-n_shift,:,i_l], v_sat[:-n_shift,:,i_l], bexp[:-n_shift,:,i_l], sathh[:-n_shift,:,i_l], hcap[:-n_shift,:,i_l], hcon[:-n_shift,:,i_l], smclsat[:-n_shift,:,i_l], tsl[:-n_shift,:,i_l], smcl[:-n_shift,:,i_l], tmax[:-n_shift,:,i_l], smclu[:-n_shift,:,i_l], smclf[:-n_shift,:,i_l], sthu[:-n_shift,:,i_l], sthf[:-n_shift,:,i_l], xice[:-n_shift,:,i_l], z_depth[:-n_shift,:,i_l], is_surface[:-n_shift,:,i_l], daily_remained_frozen[:-n_shift,:,i_l], monthly_remained_frozen[:-n_shift,:,i_l], yearly_remained_frozen[:-n_shift,:,i_l], daily_tsl[:-n_shift,:,i_l], monthly_tsl[:-n_shift,:,i_l], yearly_tsl[:-n_shift,:,i_l]
                            # assumes regular grid?

                            dx[:n_shift,:,i_l], dy[:n_shift,:,i_l], v_sat[:n_shift,:,i_l], bexp[:n_shift,:,i_l], sathh[:n_shift,:,i_l], hcap[:n_shift,:,i_l], hcon[:n_shift,:,i_l], smclsat[:n_shift,:,i_l], tsl[:n_shift,:,i_l], smcl[:n_shift,:,i_l], tmax[:n_shift,:,i_l], smclu[:n_shift,:,i_l], smclf[:n_shift,:,i_l], sthu[:n_shift,:,i_l], sthf[:n_shift,:,i_l], xice[:n_shift,:,i_l], z_depth[:n_shift,:,i_l], is_surface[:n_shift,:,i_l], daily_remained_frozen[:n_shift,:,i_l], monthly_remained_frozen[:n_shift,:,i_l], yearly_remained_frozen[:n_shift,:,i_l], daily_tsl[:n_shift,:,i_l], monthly_tsl[:n_shift,:,i_l], yearly_tsl[:n_shift,:,i_l] = dx[0,:,i_l], dy[0,:,i_l], v_sat[0,:,i_l], bexp[0,:,i_l], sathh[0,:,i_l], hcap[0,:,i_l], hcon[0,:,i_l], smclsat[0,:,i_l], tsl[0,:,i_l], smcl[0,:,i_l], tmax[0,:,i_l], smclu[0,:,i_l], smclf[0,:,i_l], sthu[0,:,i_l], sthf[0,:,i_l], xice[0,:,i_l], z_depth[0,:,i_l], is_surface[0,:,i_l], daily_remained_frozen[0,:,i_l], monthly_remained_frozen[0,:,i_l], yearly_remained_frozen[0,:,i_l], daily_tsl[0,:,i_l], monthly_tsl[0,:,i_l], yearly_tsl[0,:,i_l]

                            dx[-1,:,i_l], dy[-1,:,i_l], v_sat[-1,:,i_l], bexp[-1,:,i_l], sathh[-1,:,i_l], hcap[-1,:,i_l], hcon[-1,:,i_l], smclsat[-1,:,i_l], tsl[-1,:,i_l], smcl[-1,:,i_l], tmax[-1,:,i_l], smclu[-1,:,i_l], smclf[-1,:,i_l], sthu[-1,:,i_l], sthf[-1,:,i_l], xice[-1,:,i_l], z_depth[-1,:,i_l], is_surface[-1,:,i_l], daily_remained_frozen[-1,:,i_l], monthly_remained_frozen[-1,:,i_l], yearly_remained_frozen[-1,:,i_l], daily_tsl[-1,:,i_l], monthly_tsl[-1,:,i_l], yearly_tsl[-1,:,i_l] = dx[-2,:,i_l], dy[-2,:,i_l], v_sat[-2,:,i_l], bexp[-2,:,i_l], sathh[-2,:,i_l], hcap[-2,:,i_l], hcon[-2,:,i_l], smclsat[-2,:,i_l], tsl[-2,:,i_l], smcl[-2,:,i_l], tmax[-2,:,i_l], smclu[-2,:,i_l], smclf[-2,:,i_l], sthu[-2,:,i_l], sthf[-2,:,i_l], xice[-2,:,i_l], z_depth[-2,:,i_l], is_surface[-2,:,i_l], daily_remained_frozen[-2,:,i_l], monthly_remained_frozen[-2,:,i_l], yearly_remained_frozen[-2,:,i_l], daily_tsl[-2,:,i_l], monthly_tsl[-2,:,i_l], yearly_tsl[-2,:,i_l]

                            x_offset[i_l] += n_shift * temp_x_spacing

                for output in outputs:
                    if output['freq'] == 'daily':
                        if last_save != n_days:
                            print(f"Saving {output['freq']}_{output['vars']} month {n_months}")
                            start = time.time()
                            to_output = output['data'].isel({'time':slice(last_save,n_days)}).copy()
                            output_name = output_path + file_prefix+f"{output['freq']}_{output['vars']}_month{n_months:03}.nc"
                            output_encoding = output['encoding']
                            to_output.to_netcdf(output_name, encoding=output['encoding'])
                            end = time.time()
                            last_save = n_days
                            print(f"Saved in {end-start}s")
                        else:
                            print('Warning! last_save = n_days.') # Hopefully not needed.
                if False: # save monthly monthly
                    for output in outputs:
                        if output['freq'] == 'monthly':
                            print(f"Saving {output['freq']}_{output['vars']} month {n_months}")
                            start = time.time()
                            output['data'].to_netcdf( output_path + file_prefix+f"{output['freq']}_{output['vars']}.nc", encoding=output['encoding'])
                            end = time.time()
                            print(f"Saved in {end-start}s")
                if 'yearly' in times:
                    for output in outputs:
                        if output['freq'] in ['monthly','yearly']:
                            print(f"Saving {output['freq']}_{output['vars']} year {n_years}")
                            start = time.time()
                            output['data'].to_netcdf( output_path + file_prefix+f"{output['freq']}_{output['vars']}.nc", encoding=output['encoding'])
                            end = time.time()
                            print(f"Saved in {end-start}s")
                    n_years  += 1
                n_months +=1
            if time_to_go == 'print_daily':
                now = time.time() # Currently timing days
                elapsed = now-main_start
                if i > 0 and n_days > 0: # may result in incorrect estimate skipping day 0
                    print(f"Day {n_days} in {(now-day_start)/60:.3g} mins, {elapsed/i:.3g} s/it avg, end in {(total_timedeltas['days']-n_days)*(elapsed/n_days)/(60*60):.3g} hrs")
                day_start = time.time()
            n_days += 1
    previous_time = current_time

#----------------------
# Final save for daily
#----------------------
for output in outputs:
    if output['freq'] == 'daily':
        if last_save != n_days:
            print(f"Saving {output['freq']}_{output['vars']} final")
            start = time.time()
            output['data'].to_netcdf( output_path + file_prefix+f"{output['freq']}_{output['vars']}.nc",
                                                                            encoding=output['encoding'])
            end = time.time()
            print(f"Saved in {end-start}s")    

print('Done :)')
