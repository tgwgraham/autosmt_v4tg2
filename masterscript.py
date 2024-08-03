import os
from os import mkdir
from os.path import exists
from csbdeep.utils import normalize
import matplotlib.pyplot as plt
import skimage
from skimage import io
from skimage.measure import regionprops_table
from stardist.models import StarDist2D
import glob, shutil
import numpy as np
import time
import math
import random
import toml


macroname = 'autosmt.mac'

settings = toml.load('settings.toml')
ms = settings['macro_settings']

# Load in variables from the settings file
basefname = ms['basefname']
minarea = ms['minarea']          # minimum area of ROI
minint = ms['minint']              # minimum mean pixel intensity
maxint = ms['maxint']          # maximum mean pixel intensity
extraOC = ('extraOCs' in ms.keys())          # Whether or not to take snapshot in the extraOC channel
defaultROIsize = ms['defaultROIsize']
resizeROI = ms['resizeROI']        # Whether or not to resize the ROI to fit tightly around the cell
selectionMode = ms['selectionMode']   # options: largest - chooses the largest cell

GRIDSIZE = ms['GRIDSIZE']
STEPSIZEMICRONS = ms['STEPSIZEMICRONS']
cellFindingOC = ms['cellFindingOC']
preBleachingOC = ms['preBleachingOC']
ISOC = ms['ISOC']
bleachTime = ms['bleachTime']
ISname = ms['ISname']
INJECTAFTER = ms['INJECTAFTER']
DOINJECTION = ms['DOINJECTION'] 


MACROSTRING_0 = rf"""#define GRIDSIZE {GRIDSIZE} 
#define STEPSIZEMICRONS {STEPSIZEMICRONS} 
#define cellFindingOC "{cellFindingOC}" 
#define preBleachingOC "{preBleachingOC}"	
#define ISOC "{ISOC}"
#define bleachTime {bleachTime}
#define ISname "{ISname}"
#define INJECTAFTER {INJECTAFTER} 
int DOINJECTION = {DOINJECTION}; 

int i = 0;
int j = 0;
int counter = 0; // count of total number of FOVs (for determining injection time)
double direction = 1; 


// INITIALIZE EVERYTHING
CloseAllDocuments(0);
CameraFormatSet(1, "FMT 1x1 (X-6897) 16-bit");
CameraFormatSet(2, "FMT 1x1 (X-6897) 16-bit");
ROIEnable(0);

// LOOP OVER COLUMNS
for ( i = 0; i < GRIDSIZE ; i = i+1)
{{

// LOOP OVER ROWS
for ( j = 0; j < GRIDSIZE ; j = j+1)
{{

// INJECTION
if (DOINJECTION == 1 & counter == INJECTAFTER) {{
Int_ExecProgram("python {basefname}/syringe_inject.py");
Wait(60);
}}

// MOVE TO NEXT GRID POINT
//StgMoveXY(0.00, STEPSIZEMICRONS * direction, 1); // original version
StgMoveXY(STEPSIZEMICRONS * direction, 0.00, 1);

// TAKE INITIAL IMAGE
SelectOptConf(cellFindingOC); // Set optical configuration here for finding cells.
CameraFormatSet(1, "FMT 1x1 (X-6897) 16-bit");
CameraFormatSet(2, "FMT 1x1 (X-6897) 16-bit");
ROIEnable(0);
FineGrab();
ZoomFitToScreen();
ImageSaveAs("{basefname}/temp/snap.tif",18,0);
CloseAllDocuments(0);
"""

# Take images in other channels. These will be called by the channel names given in the settings.toml file.
MACROSTRING_1 = """// TAKE INITIAL IMAGES IN ANOTHER CHANNEL
SelectOptConf("{current_OC_name}"); // Set optical configuration here for finding cells.
CameraFormatSet(1, "FMT 1x1 (X-6897) 16-bit");
CameraFormatSet(2, "FMT 1x1 (X-6897) 16-bit");
ROIEnable(0);
FineGrab();
ZoomFitToScreen();
ImageSaveAs("{base_file_name}/temp/snap_{channel_name}.tif",18,0);
CloseAllDocuments(0);
"""

MACROSTRING_2 = rf"""// SIGNAL TO MASTERSCRIPT THAT SNAPSHOTS HAVE BEEN TAKEN
WriteFile("{basefname}/temp/done.txt","DONE",8); 

// WAIT FOR THE SNAPSHOTS TO GET PROCESSED AND MOVED BY THE PYTHON SCRIPT
while(ExistFile("{basefname}/temp/snap.tif")) {{
    "Waiting for snapshot processing.";
    }}

// MOVE TO CELL AND CONTINUE IMAGING IF A CELL IS FOUND
RunMacro("{basefname}/temp/movethere.mac");

}}

// MOVE TO NEXT ROW OF THE GRID, AND REVERSE DIRECTION
// StgMoveXY(STEPSIZEMICRONS, 0, 1); //original version
StgMoveXY(0, -STEPSIZEMICRONS, 1); 

direction = -direction;

}}
                                                                 """

def write_parent_macro(settings,basefname,macroname,extraOC=extraOC):
    with open(macroname, 'w') as fh:
        fh.write(MACROSTRING_0)
        if extraOC:
            for k in settings['macro_settings']['extraOCs']:
                fh.write(MACROSTRING_1.format(current_OC_name = settings['macro_settings']['extraOCs'][k]['extraOC'], 
                                                base_file_name = basefname, 
                                                channel_name = k))
        fh.write(MACROSTRING_2)

def movefile(basefname,snapnum,ms=ms,extraOC=extraOC):    
    print(f'Moving snap{snapnum}.')
    print('Do not close this window while the macro is running.')

    if snapnum == 1:
        currfiles = glob.glob(basefname + 'output/snaps/*')
    else:
        if snapnum==4:
            currfiles = glob.glob(basefname + 'output/snaps2/*')
        else:
            currfiles = glob.glob(basefname + 'output/snaps%d/*' % snapnum)

    maxint = 0
    for f in currfiles:
        f = f.split('\\')
        f = f[-1]
        digits = f[:-4]
        maxint = max(maxint,int(digits))

    if snapnum == 1:
        shutil.move(basefname + 'temp/snap.tif', 
                    basefname + 'output/snaps/%d.tif' % (maxint + 1))
        if extraOC:
            for k in ms['extraOCs'].keys():
                    shutil.move(basefname + f'temp/snap_{k}.tif', 
                            basefname + f'output/snaps_{k}/{maxint + 1}.tif')
        os.remove(f'{basefname}temp/done.txt')
    else:
        if snapnum==4:
            if extraOC:
                for k in ms['extraOCs'].keys():
                    shutil.move(basefname + f'temp/snap_{k}.tif', 
                            basefname + f'output/snaps2_{k}/{maxint}.tif')
            os.remove(f'{basefname}temp/done4.txt')
        else:
            shutil.move(basefname + 'temp/snap%d.tif' % snapnum, 
                        basefname + 'output/snaps%d/%d.tif' % (snapnum, maxint + 1))
            os.remove(f'{basefname}temp/done{snapnum}.txt')
        
def movesmt(basefname):
    currfiles = glob.glob(f'{basefname}/output/smt/*')
    sourcefile = glob.glob(f'{basefname}/temp/smt_temp/*nd2')

    # get last index
    maxint = 0
    for f in currfiles:
        f = f.split('\\')
        f = f[-1]
        digits = f[:-4]
        maxint = max(maxint,int(digits))

    shutil.move(sourcefile[0], f'{basefname}/output/smt/{maxint + 1}.nd2')
    os.remove(f'{basefname}temp/donesmt.txt')

def delfile(basefname,snapnum,extraOC=extraOC):
    # delete snap file rather than moving
    print(f'No cell found. Deleting snap {snapnum}.')
    print('Do not close this window while the macro is running.')

    if snapnum == 1:
        os.remove(f'{basefname}temp/snap.tif')
        os.remove(f'{basefname}temp/done.txt')
    else:
        if snapnum == 4:
            if extraOC:
                for k in ms['extraOCs'].keys():
                    os.remove(f'{basefname}temp/snap_{k}.tif')
                    if os.path.isfile(f'{basefname}/temp/done4.txt'):
                        os.remove(f'{basefname}/temp/done4.txt')
        else:
            os.remove(f'{basefname}temp/snap{snapnum}.tif')
            if os.path.isfile(f'{basefname}temp/done{snapnum}.txt'):
                os.remove(f'{basefname}temp/done{snapnum}.txt')

def locatecell(basefname,maxint,minarea,verbose=True,selectionMode=selectionMode,minint=minint,ms=ms,extraOC=extraOC):

    movestring = "StgMoveXY(%f, %f, 1);\n"
    
    # below are the pieces of the "movethere" macro
    movethere_part1 = f"""StgMoveXY(%f, %f, 1);
counter = counter + 1;
"counter";
counter;   
SelectOptConf(cellFindingOC); // Set optical configuration here for finding cells. 
FineGrab();
ZoomFitToScreen();
ImageSaveAs("{basefname}/temp/snap2.tif",18,0);
CloseAllDocuments(0);
WriteFile("{basefname}/temp/done2.txt","DONE",8);  

//WAIT FOR THE PYTHON SCRIPT TO PROCESS SNAP2
while(ExistFile("{basefname}/temp/snap2.tif")) {{
"Waiting for snap2 processing.";
}}

"""
# Take images in an extra OC
    extraOCmac = """// TAKE AN IMAGE WITH AN EXTRA OC
SelectOptConf("{current_OC}");
FineGrab();
ZoomFitToScreen();
ImageSaveAs("{base_file_name}/temp/snap_{channel_name}.tif",18,0);

"""

# Signal to masterscript that second round of snapshots have been taken in all channels
    done4mac = """WriteFile("{base_file_name}/temp/done4.txt","DONE",8);  
CloseAllDocuments(0);"""

    movethere_part2 = f"""// SHRINK ROI TO CELL AND TAKE A THIRD IMAGE
SelectOptConf(cellFindingOC); // Set optical configuration here for finding cells.
RunMacro("{basefname}/temp/roiupdate.mac");
FineGrab();
ZoomFitToScreen();
ImageSaveAs("{basefname}/temp/snap3.tif",18,0); 
WriteFile("{basefname}/temp/done3.txt","DONE",8);  
CloseAllDocuments(0);
    
// PRE-BLEACH
// Only do this if bleachTime is not equal to 0.
if (bleachTime != 0) {{
    SelectOptConf(preBleachingOC); // Set optical configuration for pre-bleaching
    CameraSet_Exposure(1, 7.000000);
    LiveSync();
    Wait(bleachTime);    
    Freeze();
}}

// RUN ILLUMINATION SEQUENCE
SelectOptConf(ISOC); // Set optical configuration for illumination sequence
CameraSet_Exposure(1, 7.000000);
IlluminationSequence_Run(ISname);    
CloseAllDocuments(0);
WriteFile("{basefname}/temp/donesmt.txt","DONE",8);  

// MOVE BACK TO GRID POINT CENTER
StgMoveXY(%f, %f, 1);"""

    # creates a pretrained model
    model = StarDist2D.from_pretrained('2D_versatile_fluo')

    #I = plt.imread('../output/snaps/3.tif')
    I = io.imread(f"{basefname}/temp/snap.tif", plugin='pil')
    labels, details = model.predict_instances(normalize(I), prob_thresh=0.5) 

    regiondata = skimage.measure.regionprops(labels)

    # get a bunch of properties of the ROIs
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

    background = I[labels==0].mean()
    
    sel = ((rp['intensity_mean']-background<maxint) & (rp['area'] > minarea) & (rp['intensity_mean']-background>minint))

    # If additional snapshot images were taken in other OCs, then get region properties and gate on intensity in those channels
    if extraOC:
        rp_extra = {}
        back_extra = {}
        eos = ms['extraOCs']
        for k in eos.keys():
            minint2 = eos[k]['minint']
            maxint2 = eos[k]['maxint']
            I2 = io.imread(f"{basefname}/temp/snap_{k}.tif", plugin='pil')
            labels2, details2 = model.predict_instances(normalize(I2), prob_thresh=0.5) 
            rp_extra[k] = regionprops_table(labels, I2,
                      properties=('intensity_mean',)
                          )
            back_extra[k] = I2[labels2==0].mean()
            # select only the subset of cells that fit within the bounds
            sel = sel & (rp_extra[k]['intensity_mean']-back_extra[k]<maxint2) & (rp_extra[k]['intensity_mean']-back_extra[k]>minint2)
    
    # print realtime intensities of cells, etc.
    if verbose:
        if extraOC:
            print('intensity_mean\t\tarea',end='')
            for k in eos.keys():
                print(f"\t\tintensity_mean_{k}",end='')
            print()
            for i in range(len(regiondata)):
                print(f"{rp['intensity_mean'][i]-background}\t\t{rp['area'][i]} ",end="")
                for k in eos.keys():
                    print(f"\t\t{rp_extra[k]['intensity_mean'][i]-back_extra[k]}",end='')
                if sel[i]:
                    print("***")
                else:
                    print()
        else:
            print('intensity_mean\t\tarea')
            for i in range(len(regiondata)):
                print(f"{rp['intensity_mean'][i]-background}\t\t{rp['area'][i]} ",end="")
                if sel[i]:
                    print("***")
                else:
                    print()

    large_regions = [];
    for i,region in enumerate(regiondata):    
      if sel[i]:
        large_regions.append(region)
        
    # if verbose:
        # print(large_regions)

    if not bool(large_regions):
      shiftx = 0
      shifty = 0
      with open(f'{basefname}/temp/movethere.mac', 'w') as fh:
        fh.write(movestring % (shiftx,-shifty))     
      return False # return value indicating that a suitable cell was not found

    else:

      # "largest_circular_region" is a bit of a misnomer, because I've added other possible selection modes
      if selectionMode == "largest":
        largest_circular_region = max(large_regions, key=lambda item: item.area)
      elif selectionMode == "random":
        largest_circular_region = max(large_regions, key=lambda item: random.random())
      else: # defaults to largest nucleus
        print("selectionMode not recognized! Defaulting to largest ROI.")
        largest_circular_region = max(large_regions, key=lambda item: item.area)

      centroid = largest_circular_region.centroid

      # remember that we used "valid" mode in the convolution, so indices are already shifted by half of the box width
      micronsperpx = 0.16 
      xcenter = centroid[1]
      ycenter = centroid[0]
      shiftx = micronsperpx * (xcenter - 256)
      shifty = micronsperpx * (ycenter - 256)
      
      
      # New version: Use the strings above
      with open(f'{basefname}/temp/movethere.mac', 'w') as fh:
        if extraOC:
            fh.writelines(movethere_part1 % (shiftx,shifty)) 
            for k in ms['extraOCs'].keys():
                fh.writelines(extraOCmac.format(current_OC=ms['extraOCs'][k]['extraOC'],channel_name=k,base_file_name=basefname)) 
            fh.writelines(done4mac.format(base_file_name=basefname))
            fh.writelines(movethere_part2 % (-shiftx,-shifty))
        else:
            fh.writelines(movethere_part1 % (shiftx,shifty) + movethere_part2 % (-shiftx,-shifty))
      return True # return value indicating that a suitable cell was found


def relocatecell(basefname,verbose=True,
        resizeROI=resizeROI,selectionMode=selectionMode,defaultROIsize=defaultROIsize):
    # creates a pretrained model
    model = StarDist2D.from_pretrained('2D_versatile_fluo')

    I = io.imread(f"{basefname}/temp/snap2.tif", plugin='pil')
    labels, details = model.predict_instances(normalize(I), prob_thresh=0.5) 

    # get a bunch of properties of the ROIs
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

    background = I[labels==0].mean()


    sel = ((rp['centroid-0'] > 199)
                & (rp['centroid-0'] < 312)
                & (rp['centroid-1'] > 199)
                & (rp['centroid-1'] < 312)
                & (rp['area'] > minarea))

    regiondata = skimage.measure.regionprops(labels)

    large_regions = []
    for i,region in enumerate(regiondata):    
      if sel[i]:
        large_regions.append(region)    
    
    if verbose:
        print(large_regions)

    if bool(large_regions):

        # "largest_circular_region" is a bit of a misnomer, because I've added other possible selection modes
        if selectionMode == "largest":
            largest_circular_region = max(large_regions, key=lambda item: item.area)
        elif selectionMode == "random":
            largest_circular_region = max(large_regions, key=lambda item: random.random())
        else: # defaults to largest nucleus
            print("selectionMode not recognized! Defaulting to largest ROI.")
            largest_circular_region = max(large_regions, key=lambda item: item.area)

        bbox = largest_circular_region.bbox
        centroid = largest_circular_region.centroid

        if resizeROI:
            minx = bbox[1]
            maxx = bbox[3]
            miny = bbox[0]
            maxy = bbox[2]
        else:
            meanx = (bbox[1] + bbox[3])/2
            meany = (bbox[0] + bbox[2])/2
            minx = meanx - math.floor(defaultROIsize[0]/2)
            maxx = meanx + math.ceil(defaultROIsize[0]/2)
            miny = meany - math.floor(defaultROIsize[1]/2)
            maxy = meany + math.ceil(defaultROIsize[1]/2)

    else:  # if a cell is not relocated, set the ROI by default to be a default sized ROI in the center of the FOV
        minx = 256 - math.floor(defaultROIsize[0]/2)
        maxx = 256 + math.ceil(defaultROIsize[0]/2)
        miny = 256 - math.floor(defaultROIsize[1]/2)
        maxy = 256 + math.ceil(defaultROIsize[1]/2)
        
        if verbose:
            print('Cell not relocated. Defaulting to centered 150 x 150 px ROI.')
        
    ROImacrostring = """ROISet(%d,%d,%d,%d);
    CameraFormatSet(1, "FMT 1x1 (X-7244) 16-bit");
    CameraFormatSet(2, "FMT 1x1 (X-7244) 16-bit");
    ROIEnable(1); """ % (minx,miny,maxx,maxy) 

    with open(f'{basefname}/output/rois.txt', 'a') as fh:
        fh.write('%d,%d,%d,%d\n' % (minx,miny,maxx,maxy))

    with open(f'{basefname}/temp/roiupdate.mac', 'w') as fh:
        fh.write(ROImacrostring)


if __name__ == "__main__":

    write_parent_macro(settings,basefname,macroname)
    
    
    # make output directories
    outdirs = [f'{basefname}/{j}/' for j in ['output',
        'output/snaps',
        'output/snaps2',
        'output/snaps3',
        'output/smt',
        'output/scripts',
        'temp',
        'temp/smt_temp']]

    for outdir in outdirs:
        if not exists(outdir):
            mkdir(outdir)
        
    # make snaps folders for additional OCs
    if extraOC:
        for k in ms['extraOCs'].keys():
            outdir = f'{basefname}/output/snaps_{k}/'
            if not exists(outdir):
                mkdir(outdir)
            outdir = f'{basefname}/output/snaps2_{k}/'
            if not exists(outdir):
                mkdir(outdir)
    
    # copy all scripts, including this one, to the output/scripts folder for future reference
    filestocopy = []
    for ext in ['py','toml','rtf','mac','ipynb']:
        filestocopy.extend(glob.glob(f'*.{ext}'))
    print(filestocopy)
    for f in filestocopy:
        shutil.copy(f,'output/scripts/')
        
    # loop that will keep running to detect, process, and move files as they are created by NIS Elements 
    while True:
        if os.path.isfile(f'{basefname}/temp/done.txt'):
            # attempt to find a cell in the first snapshot. If it is found, move the snapshot and write out movethere.mac. 
            # Otherwise, delete the snapshot and write out a placeholder movethere.mac
            if locatecell(basefname,maxint,minarea):
                movefile(basefname,1)
            else:
                delfile(basefname,1)
                delfile(basefname,4)
        if os.path.isfile(f'{basefname}/temp/done2.txt'):
            # relocate cell in second snapshot and move second snapshot to snaps2 folder
            relocatecell(basefname)
            movefile(basefname,2)
        if os.path.isfile(f'{basefname}/temp/done3.txt'):
            try:
                movefile(basefname,3)
            except Exception as e:
                print(e)
        if os.path.isfile(f'{basefname}/temp/done4.txt'):
            try:
                movefile(basefname,4)
            except Exception as e:
                print(e)
        if os.path.isfile(f'{basefname}/temp/donesmt.txt'):
            movesmt(basefname)
        time.sleep(0.1)
    
# TODO: 
# Add a stop signal from the macro to break out of the while loop, clean up temp files, and exit python program














