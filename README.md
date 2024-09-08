# White_2024
Scripts for analysis and figure generation in White et al., 2024. 

## Data download
The single molecule data contained in White et al., 2024 (raw video files and simplified .mat data structures) can be downloaded from:  Note, all single molecule videos were collected using MicroManager 2.0 

## For White 2024:
1. Download/clone this repo. In Matlab, add this repo to your file path.
2. To analyze data from White et al., 2024, additionally download the data linked above and add the .mat to file path. Many scripts in "plot_figures" were hardcoded to the exact file path on my Mac, so will need to edit this to run these scripts.

## General Features: 
The GUIs contained in this repo (smtoolbox2024) were written to analyze 1 or 2 color co-localization single molecule spectroscopy (aka CoSMoS) vidoes.
* 'smVideoProcessing' loads in .tif stacks from MicroManager 2.0 and can be used to align images, correct drift, detect single molecules, and intergerate single molecules across frames. The output is a .mat file that can be loaded by 'smTraceViewer'
* 'smTraceViewer' is for the visulization and idealization if (multi-state) single molecule data. The output is an upadted .mat data structure that is then loaded directly into Matlab for custom analyis (see plot_figures for some examples)
* smtoolbox_2024/utilties contains a variety of useful (IMO) functinos for single molecule analysis. For example, 'fitDwells' can fit dwell times to single and biexponetial distributions using maximum likelihood estimation and 'plotDwells' is a quick way to visualize a histogram of dwell times overlaid with the MLE result.  

## Concluding Remarks
This repo exists to upload my local scripts and codebase used in White 2024. The code contained herein was not intended for distribution in this state, but only for clairty in how data were treated and how figures were made. As life would have it, I have moved on to industry as of Septemeber 2023 and am unlikely to update, document, or respond to issues in this repo. But who knows- maybe I will get bored or venture back one day. 




