# ADV Data Processing in MATLAB

There is currently one high level script in the repository for ADV data post-processing (adv_post_process.m). The script reads ASCII file(s) from a Nortek Vectrino and calls multiple functions to post-process the ADV data.

OBJECTIVES:
- Arrange the data into appropriate arrays
- Apply manual data windowing based on 2D traverse motion parameters
- Determine inferred pitch and yaw of probe from 3D velocity measurements
- Apply despiking algorithm of Goring and Nikora (2002)
- Apply rotational corrections to the data to account for any probe misalignment (code apdapted from Mujal, et al. UPC (2015) https://github.com/vicente-medina/TurbulenceToolbox)
- Determine noise floor in velocity signal
- Save data


INSTRUCTIONS:
- The user is to add their own filepaths in all .m files, otherwise errors will result
- Press Run and follow instructions in the Command Window


NOTES:
- Vertical and cross-stream measurement profiles contain a measurement increment of 100 mm
- Vertical profiles contain nine measurement points throughout the 1000 mm channel depth
- Horitzontal profiles contain ten measurement points across the 1100 mm channel width
- Transect measurement increment is 125 mm in both the horizontal and vertical directions
- Measurements performed on a 2D traverse, the speed and acceleration parameters for which can be modified in profiling_post_processing.m


TO FOLLOW:
- A script for ADV data analysis and its associated functions will follow at a later date
