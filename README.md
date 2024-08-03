# autosmt_v4tg2

Code for single-molecule tracking (SMT) experiment automation. The Tjian-Darzacq group operates this on a Nikon Ti-E microscope running NIS Elements 4.60. Here is a overall description of how this works:

1. Edit settings.toml to define your experiment. Details of what each variable does are given in the comments. Snapshots can be taken in an arbitrary number of channels by defining additional
macro_settings.extraOCs entries. The microscope will raster over a GRIDSIZE x GRIDSIZE grid with spacing of STEPSIZEMICRONS. Be careful that GRIDSIZE * STEPSIZEMICRONS does not exceed the size of your imaging dish in microns. 
2. Run masterscript.py in a terminal, and keep it running throughout the experiment. This will write out a macro in the same folder called "v4tg2.mac". (You can change this name as desired by changing the "macroname" variable at the top of masterscript.py)
3. Optional: Run realtime_analysis.py in a separate terminal, and keep it running throughout the experiment. realtime_analysis.py tracks molecules in real time using quot, which can be found here: https://github.com/alecheckert/quot. It also makes plots of the number of localizations per frame, the spatial distribution of localizations, and the step size distribution. Running this script is not essential, however.
4. Define your illumination sequence in the Illumination Sequence panel of NIS Elements, and be sure to set it to save to the folder "{basefname}/temp/smt_temp/", where you should replace {basefname} with whatever base folder name you listed in settings.toml. Launch v4tg2.mac from NIS elements.
5. Double-check the optical configurations in NIS Elements that you have specified in the "cellFindingOC", "preBleachingOC", "ISOC", and "extraOC" variables of settings.toml.
6. Focus on your cells in NIS Elements, and launch macro v4tg2.mac.

You can pause and resume the macro as needed to refocus the cells, manually add compounds to the sample, etc. If you have to abort the macro, you may have to do this twice, as it runs a few nested macros-within-a-macro. It is usually best to re-start masterscript.py and realtime_analysis.py.

The data and analysis will be in a folder called "output". It is easiest to transfer this after the experiment to another folder on the same hard drive and rename it to something more informative. You can delete the "temp" folder after the experiment is over.

Notebook realtime_PAPA_stats.ipynb can be used to visualize PAPA results in real time.
