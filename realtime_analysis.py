from quot import read_config, track_file
import glob
from os.path import exists
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os
import time
from csbdeep.utils import normalize
from skimage import io
from stardist.models import StarDist2D
import re
from skimage.measure import regionprops_table
import imageio
import toml

config = toml.load('settings.toml')
ms = config['macro_settings']

basefname = ms['basefname']
#config = read_config(f'{basefname}settings.toml')
nd2folder = f'{basefname}/output/smt'
# extra_snap_folders = [] # list of folders containing extra snapshots
# ocnames = []            # names of extra snapshot channels
if 'extraOCs' in ms.keys():
    ocnames = list(ms['extraOCs'].keys())
    extra_snap_folders = [f"{basefname}/output/snaps2_{k}/" for k in ocnames]
else:
    ocnames = []
    extra_snap_folders = []

ispapa = ms['ispapa']

umperpx = 0.16


def get_gv_framenum(config):
    """get_gv_framenum(config)
    
    returns arrays containing the frame indices of frames before and after green and violet excitation
    (PAPA and DR, respectively)"""
    
    def add_fw_and_flatten(x,ncycles,fwrange):
        rv = x.reshape(ncycles,1) + fwrange
        rv = rv.flatten()
        return rv
    
    IS = config['illumination_sequence']
    ncycles = IS['ncycles']
    r = IS['r']
    v = IS['v']
    g = IS['g']
    fw = IS['framewindow']
    gfirst = IS['gfirst']
    
    t = 4*r + v + g
    
    fwrange = np.arange(0,fw)
    
    if gfirst: # if the green pulse occurs first in the cycle
        gpre = np.arange(0,ncycles)*t # indices of frames in frame window before green
        gpost = np.arange(0,ncycles)*t + r + g # indices of frames in frame window after green
        vpre = np.arange(0,ncycles)*t + 2*r + g # indices of frames in frame window before violet
        vpost = np.arange(0,ncycles)*t + 3*r + v + g # indices of frames in frame window after violet
    else: # if the violet pulse occurs first in the cycle
        vpre = np.arange(0,ncycles)*t # indices of frames in frame window before violet
        vpost = np.arange(0,ncycles)*t + r + v # indices of frames in frame window after violet
        gpre = np.arange(0,ncycles)*t + 2*r + v # indices of frames in frame window before green
        gpost = np.arange(0,ncycles)*t + 3*r + v + g # indices of frames in frame window after green
    
    gpre = add_fw_and_flatten(gpre,ncycles,fwrange)
    gpost = add_fw_and_flatten(gpost,ncycles,fwrange)
    vpre = add_fw_and_flatten(vpre,ncycles,fwrange)
    vpost = add_fw_and_flatten(vpost,ncycles,fwrange)
    
    return({'vpre':vpre,'vpost':vpost,'gpre':gpre,'gpost':gpost})


gvfn = get_gv_framenum(config)
gframes = gvfn['gpost']
vframes = gvfn['vpost']

# # Post-green and post-violet frame numbers for 30-frame PAPA protocol
# gframes = np.array([])
# for j in range(5):
    # gframes = np.append(gframes, np.arange(90,120) + 120*j)
# vframes = np.array([])
# for j in range(5):
    # vframes = np.append(vframes, np.arange(30,60) + 120*j)

columns = ['y','x','I0','bg','y_err','x_err','I0_err','bg_err','H_det','error_flag', 
      'snr','rmse','n_iter','y_detect','x_detect','frame','loc_idx','trajectory',  
      'subproblem_n_traj','subproblem_n_locs']

os.makedirs(f'{basefname}/output/tracking',exist_ok=True)
os.makedirs(f'{basefname}/output/masks',exist_ok=True)
os.makedirs(f'{basefname}/output/roi_measurements',exist_ok=True)
os.makedirs(f'{basefname}/output/displacements',exist_ok=True)

model = StarDist2D.from_pretrained('2D_versatile_fluo')


def get_ND2_time(fname):
        
    """
    Get timestamp in Julian Date Number from a Nikon nd2 file.
    
    args
        fname  :   string - path to the nd2 file
        
    returns
        Julian Date Number as a floating point number 
    
    """
    
    bytesfromend = int(1e6) # where to start reading in the metadata relative the end of the file
    
    # Julian Date Number pattern to match in the metadata at the end of the ND2 file
    pattern = r"<ModifiedAtJDN runtype=\"double\" value=\"(.+?)\"/>" 

    with open(fname, "rb") as file:
        # Read the bytes at the end of the file
        file.seek(os.path.getsize(fname)-bytesfromend)
        bytes = file.read()
        decoded = bytes.decode(errors='ignore')
        numbers = re.findall(pattern, decoded)
        if not numbers:
            return float('nan')
        else:
            return float(numbers[0])


def getalldisp(trajs):
    """
getalldisp(trajs)

returns a list of all single-molecule displacements from a csv file containing
single-molecule trajectories. Accepts a pandas dataframe containing trajectory data

This function assumes that the x and y coordinates are in the first two columns
of the csv file.
    """
    trajn = trajs['trajectory'].unique()
    sz = 0
    for j in range(trajn.size):
        tcount = np.sum(trajs['trajectory'] == trajn[j])
        if tcount > 1:
            sz = sz + tcount - 1
    rv = np.zeros(sz)
    counter = 0
    for j in range(trajn.size):
        selector = (trajs['trajectory'] == trajn[j])
        if np.sum(selector) > 1:
            currtraj = trajs[selector][['x','y']].to_numpy()
            rv[counter:(counter + currtraj.shape[0] - 1)] = \
                np.sqrt(np.power(currtraj[1:, 0] - currtraj[:-1, 0], 2) 
                        + np.power(currtraj[1:, 1] - currtraj[:-1, 1], 2))
            counter = counter + currtraj.shape[0] - 1
    return rv

# The loop below runs continuously. The the first step, it locates snapshot files to be segmented and 
# measured and does this first, because it is faster. In the second step, it locates ND2 files to 
# be tracked and tracks them.

while True:
    # step 1: get a list of all snapshot files that have not yet been segmented and measured
    flist = glob.glob(f"{basefname}/output/snaps2/*.tif")
    # include files in the file list for which snapshots in all channels are available
    for folder in extra_snap_folders:
        flist = [f for f in flist if exists(f"{folder}/{os.path.basename(f)}")]
    flist = [f for f in flist if not exists(f"{basefname}/output/masks/{os.path.basename(f)[:-4]}.csv")]
    for f in flist:

        # output of image segmentation
        segfname = f'{basefname}/output/masks/{os.path.basename(f)[:-4]}.csv'
        # output CSV containing ROI properties
        propfname = f'{basefname}/output/roi_measurements/{os.path.basename(f)[:-4]}.csv'
        
        
        I = io.imread(f, plugin='pil')
        labels, details = model.predict_instances(normalize(I), prob_thresh=0.5) 
        np.savetxt(segfname,labels,delimiter=',',fmt='%d')
        
        # get region properties of all the ROIs in the image
        rp = regionprops_table(labels, I,
                      properties=('centroid',
                                 'orientation',
                                 'axis_major_length',
                                 'axis_minor_length',
                                  'area',
                                  'area_filled',
                                  'eccentricity',
                                  'intensity_mean',
                                 )
                          )
        df = pd.DataFrame(rp,index = range(1,len(rp['area'])+1))
        df['background'] = I[labels==0].mean()

        # measure intensity mean and background in each snapshot channel, and append this as another column of the dataframe
        # Note that "background" for the other channels may be meaningless if there is signal in the other channels outside
        # of the nuclear masks defined by the first channel
        for i,folder in enumerate(extra_snap_folders):
            # get intensity mean for each addition channel
            
            # read in image for this channel                
            I = io.imread(f"{folder}/{os.path.basename(f)}", plugin='pil')
            # measure only mean intensity
            rp = regionprops_table(labels, I,properties=('intensity_mean',))
            newcolumns = pd.DataFrame(rp,index = range(1,len(rp['intensity_mean'])+1))
            # get "background" intensity outside of all ROIs. Note that this may not be very meaningful
            # unless the segmentation in the first channel is valid for all channels
            newcolumns[f'background_{ocnames[i]}'] = I[labels==0].mean()
            # the intensity in each channel will be called by whatever name it was given in the macro settings file
            newcolumns.rename(columns={'intensity_mean':f'intensity_mean_{ocnames[i]}'},inplace=True)
            df = pd.concat([df,newcolumns],axis=1)
            
        # get timestamp for nd2 file associated with this snapshot file
        nd2fname = f'{nd2folder}/{os.path.basename(f)[:-4]}.nd2'
        
        # try getting timestamp until we succeed (i.e., file exists and is available for reading)
        gottime = False
        while not gottime:
            try:
                df['time'] = get_ND2_time(nd2fname)
                gottime = True
            except:
                pass
        
        
        # Write out properties to a csv file.
        df.to_csv(propfname)
        
    # step 2: get a list of all ND2 files that have not yet been tracked
    flist = glob.glob(f"{basefname}/output/smt/*.nd2")
    flist = [j for j in flist if not exists(f'{basefname}/output/tracking/{os.path.split(j)[-1][:-4]}.csv')]

    # loop over all ND2 files
    for f in flist:
        trajs = track_file(f, **config)
        f2 = os.path.split(f)
        f2 = f2[-1][:-4]
        try:
            trajs.to_csv(f'{basefname}/output/tracking/{f2}.csv', index=False, columns=columns)
        except:
            with open(f'{basefname}/output/tracking/{f2}.csv','w') as fh:
                [fh.write('f{j},') for j in columns]
            continue

        if ispapa:
            # output files containing green and violet displacements
            gdispfname = f'{basefname}/output/displacements/{os.path.basename(f)[:-4]}_gdisp.csv'
            vdispfname = f'{basefname}/output/displacements/{os.path.basename(f)[:-4]}_vdisp.csv'
        else:
            dispfname = f'{basefname}/output/displacements/{os.path.basename(f)[:-4]}_disp.csv'


        fig,ax = plt.subplots(2,2,figsize=(10,6))
        maxframe = trajs['frame'].max()
        counts, bins = np.histogram(trajs['frame'].to_numpy(),bins=np.arange(0,maxframe+1)-0.5)
        heatmap, xedges, yedges = np.histogram2d(
                np.arange(0,maxframe), counts, 
                bins=[int(maxframe/5),np.arange(0,max(counts)+1,1)])
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        #plt.figure(figsize=[6,6])
        ax[0,0].imshow(heatmap.T, extent=extent, origin='lower',aspect=maxframe/max(counts),cmap='gray_r',interpolation='none')
        ax[0,0].set_xlabel('Frame number')
        ax[0,0].set_ylabel('Number of localizations')

        im = imageio.imread(f'{basefname}/output/snaps3/{f2}.tif')
        
        if ispapa:
            ax[0,1].imshow(im,cmap='gray')

            sel = trajs['frame'].isin(gframes)
            currtrajs = trajs[sel]
            for tnum in currtrajs['trajectory'].unique():
                currtraj = currtrajs[currtrajs['trajectory']==tnum]
                ax[0,1].plot(currtraj['x'],currtraj['y'],'-',color='green')
                ax[0,1].plot(currtraj['x'].iloc[0],currtraj['y'].iloc[0],'.',color='darkgreen')
                
            sel = trajs['frame'].isin(vframes)
            currtrajs = trajs[sel]
            for tnum in currtrajs['trajectory'].unique():
                currtraj = currtrajs[currtrajs['trajectory']==tnum]
                ax[0,1].plot(currtraj['x'],currtraj['y'],'-',color='violet')
                ax[0,1].plot(currtraj['x'].iloc[0],currtraj['y'].iloc[0],'.',color='darkviolet')
            ax[0,1].axis('off');
            
            
        
            gtraj = trajs[trajs['frame'].isin(gframes)]
            gd = getalldisp(gtraj)
            vtraj = trajs[trajs['frame'].isin(vframes)]
            vd = getalldisp(vtraj)
            gh,gx = np.histogram(gd*umperpx,bins=50);
            gh = gh/gh.sum()
            vh,vx = np.histogram(vd*umperpx,bins=50);
            vh = vh/vh.sum()
            
            np.savetxt(gdispfname,gd,delimiter=',')
            np.savetxt(vdispfname,vd,delimiter=',')

            ax[1,0].plot((gx[1:]+gx[:-1])/2,gh,color='darkgreen')
            ax[1,0].plot((vx[1:]+vx[:-1])/2,vh,color='darkviolet')
            ax[1,0].set_xlabel('Displacement (µm)')
            ax[1,0].set_ylabel('Frequency');
            
        
        else:
            ax[0,1].imshow(im,cmap='gray')
            ax[0,1].plot(trajs['x'],trajs['y'],'.',color='blue')

            d = getalldisp(trajs)
            gh,gx = np.histogram(d*umperpx,bins=50);
            gh = gh/gh.sum()
            ax[1,0].plot((gx[1:]+gx[:-1])/2,gh,color='blue')
            ax[1,0].set_xlabel('Displacement (µm)')
            ax[1,0].set_ylabel('Frequency');
            
            np.savetxt(dispfname,d,delimiter=',')
            
        ax[1,1].set_visible(False)
        
        fig.savefig(rf'{basefname}/output/tracking/{f2}.png')
        plt.close()
    time.sleep(0.1)