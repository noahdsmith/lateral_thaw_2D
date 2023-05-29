import xarray as xa
import numpy as np
from tqdm import tqdm
"""
Plotting
"""
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature
import matplotlib as mpl
import geopandas as gpd
import datetime
import calendar
from scipy import interpolate

# Generate colourmap
l = 128
from matplotlib import cm
from matplotlib.colors import ListedColormap
hot2 = cm.get_cmap('jet')(np.linspace(0.64,1,l))
cold2 = cm.get_cmap('jet')(np.linspace(0,0.4,l))
t_arctic2 = ListedColormap(np.vstack((cold2,hot2)), name='t_arctic2')

degree_sign= u'\N{DEGREE SIGN}'
temp_label = 'Temperature ('+ degree_sign +'C)'

from matplotlib.animation import FuncAnimation
from IPython import display
from IPython.display import Video

# Generate colourmap
l = 128
a, b, c = np.linspace(1,0,l), np.ones(l), np.zeros(l)

sub_zero = ListedColormap(np.dstack((c,a,b,b))[0], name = 'sub_zero')

hot = cm.get_cmap('autumn')(np.linspace(1,0,l))
cold = sub_zero(np.linspace(1,0,l))

t_arctic = ListedColormap(np.vstack((cold,hot)), name='t_arctic')

hot2 = cm.get_cmap('jet')(np.linspace(0.64,1,l))
cold2 = cm.get_cmap('jet')(np.linspace(0,0.4,l))

t_arctic2 = ListedColormap(np.vstack((cold2,hot2)), name='t_arctic2')

ice2 = ListedColormap(cold2, name='ice2')

purple_pink_white = cm.get_cmap('RdPu_r')(np.linspace(1,0,l))
hot_inv = cm.get_cmap('autumn')(np.linspace(0,1,l))

not_blue = ListedColormap(np.vstack((purple_pink_white,hot_inv)), name='not_blue')

hot_crop = cm.get_cmap('gist_heat')(np.linspace(0.8,0,l))
magma_crop = cm.get_cmap('gnuplot2')(np.linspace(0,0.6,l))                    
hotblackcool = ListedColormap(np.vstack((hot_crop,magma_crop)), name='hotblackcool')

def plot_and_save(filename = 'test_output/10_global_offset__69.25N_25.25Emonthly_full.nc', 
                  outputname = 'test_output_xice_interpolate_sthu2_without_nanfilter.mp4'):
    test_out = xa.open_mfdataset(filename)
    bob = test_out
    threshold = 0.15
    bob['thaw_positions'] = (['time'],np.array([interpolate.interp1d(bob['z_depth_surface'][nt,::-1,0], bob.x[::-1], bounds_error=False, assume_sorted = True)(threshold) for nt in range(len(bob.time))]))
    plt.ioff()
    maxlayer = 20
    #new!
    fig = plt.figure(constrained_layout=True,figsize = (15,6))
    gs = fig.add_gridspec(3, 4)
    ax1 = fig.add_subplot(gs[:,:2])
    ax2 = fig.add_subplot(gs[0, 2])
    ax3 = fig.add_subplot(gs[1, 2])
    ax4 = fig.add_subplot(gs[2, 2])
    ax5 = fig.add_subplot(gs[0, 3])
    ax6 = fig.add_subplot(gs[1:, 3])
    
    T_cropped = bob.t_soil[:,:,:maxlayer,0]
    vmin,vmax = np.nanmin(T_cropped.data).compute(), np.nanmax(T_cropped.data).compute()
    norm = mpl.colors.TwoSlopeNorm(vmin = vmin,vcenter=0,vmax=vmax)
    im1 = T_cropped[0,:,:].plot(x='x',ax=ax1,norm=norm,cmap=t_arctic2)
    surface_line = ax1.plot(bob.x,bob.z_depth_surface[0,:,0],color='green')
    #im1 = ax1.imshow(np.transpose(tsl_list[0]),cmap=t_arctic2,norm=norm)
    #fig.colorbar(im1,ax=ax1)
    ax1.set_title('t_soil')
    
    xice_cropped = bob.xice[:,:,:maxlayer,0]
    vmax = np.nanmax(xice_cropped.data).compute()
    im2 = xice_cropped[0,:,:].plot(x='x',ax=ax2,vmin=0,vmax=vmax,cmap='inferno')
    ax2.set_title('xice')
    
    sthuf_cropped = bob.sthuf[:,:,:maxlayer,0]
    im3 = sthuf_cropped[0,:,:].plot(x='x',ax=ax3,vmin=0,vmax=1,cmap='Blues')
    ax3.set_title('sthuf')
    
    sthf_cropped = bob.sthf[:,:,:maxlayer,0]
    im4 = sthf_cropped[0,:,:].plot(x='x',ax=ax4,vmin=0,vmax=1,cmap='Blues')
    ax4.set_title('sthf')
    
    surface_line_5 = ax5.plot(bob.x,bob.z_depth_surface[0,:,0],color='green',label='surface')
    xice_line = ax5.plot(bob.x,(bob.xice[0,:,:,0]*bob.dy[0,:,:,0]).sum(axis = 1),color='grey',label='xice total')
    ax5.legend()
    ax5.set_title('surface and xice total')

    thaw_line_6 = ax6.plot(bob.time.data,bob.thaw_positions.data)
    ax6.set_title('{threshold}m elevation position')
    
    fps = 5
    skip=1
    def animate(i):
        im1.set_array(np.transpose(T_cropped[i*skip,:,:].data))
        surface_line[0].set_data(bob.x,bob.z_depth_surface[i*skip,:,0].data)
        im2.set_array(np.transpose(xice_cropped[i*skip,:,:].data))
        im3.set_array(np.transpose(sthuf_cropped[i*skip,:,:].data))
        im4.set_array(np.transpose(sthf_cropped[i*skip,:,:].data))
        surface_line_5[0].set_data(bob.x,bob.z_depth_surface[i*skip,:,0].data)
        xice_line[0].set_data(bob.x,(bob.xice[i*skip,:,:,0]*bob.dy[i*skip,:,:,0]).sum(axis = 1).data)
        thaw_line_6[0].set_data(bob.time[:i*skip].data,bob.thaw_positions[:i*skip].data)
        time = T_cropped.time[i*skip]
        ax1.set_title(f"Frame: {i*skip}")
        ax2.set_title(f"Time: {time.dt.date.data})")
    #fig.tight_layout()
    anim = FuncAnimation(fig, animate, interval = 1000/fps, frames = int(len(T_cropped.time.data)-12/skip)-1)
    video=anim.to_html5_video()
    html = display.HTML(video)
    display.display(html)
    plt.close()
    anim.save(outputname,writer=mpl.animation.FFMpegWriter(fps=10))
    bob.close()
    Video(outputname)

from global_lat_thaw_parm_mod import output_path, prefix_prefix, sp

# print(f"Plotting {output_path}{prefix_prefix}{sp.prefix_suffix}monthly_full.nc")
# print(f"Saving as {output_path}{prefix_prefix}{sp.prefix_suffix}_video.mp4")

# plot_and_save(filename = f"{output_path}{prefix_prefix}{sp.prefix_suffix}monthly_full.nc", 
#                   outputname = f"{output_path}{prefix_prefix}{sp.prefix_suffix}_video.mp4")

plot_and_save(filename = f"{output_path}{prefix_prefix}_{69.25}N_{25.25}E_{sp.param_option}monthly_full.nc", 
                  outputname = f"{output_path}{prefix_prefix}_{69.25}N_{25.25}E_{sp.param_option}_video.mp4")
