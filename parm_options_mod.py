parm_options_dict = {'default':{}} # The default one is blank, so it doesn't overwrite any of the default parameters.
i = 0
for sel_xice in [0.1,0.3,0.55]: # 0.1,0.25,0.4
    for sel_mineral_depth in [0.75,1.0,1.5]:# 1.5, 0.75
        for sel_elevation in [0.7, 2,3]: # 0.6, 2
            for sel_hist_or_sspno in ['historical', 'ssp370']:# 'historical', 'ssp126' , 'ssp370', 'ssp585'
                parm_options_dict[i] = dict(
                hist_or_sspno = sel_hist_or_sspno, # for future runs [ ssp585, ssp370, ssp126]
                mineral_depth_bs = sel_mineral_depth,
                xice_content = sel_xice,
                xice_depth_bs = sel_mineral_depth, # Usually set to top of mineral layer.
                high_bit_elevation = sel_elevation,
                param_option = i,
        		depth_beyond_7 = 14,
        		grid_option = 'dense_longer'
                )
                i += 1
