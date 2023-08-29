# ADV Data Processing using MATLAB

There is currently one high level script in the repository for ADV data post-processing (adv_post_process.m). The script reads ASCII file(s) from a Nortek Vectrino and calls multiple functions to post-process the ADV data.

OBJECTIVES:
- Arrange the data into ordered arrays
- Apply manual data windowing based on 2D traverse motion parameters
- Determine inferred pitch and yaw of probe from 3D velocity measurements
- Apply despiking algorithm of Goring and Nikora (2002) with Wahl (2003) amendments
- Apply rotational corrections to the data to account for any probe misalignment (code adapted from Mujal, et al. UPC (2015) https://github.com/vicente-medina/TurbulenceToolbox)
- Determine noise floor in velocity signal
- Save data


INSTRUCTIONS:
- The user is to add their own filepaths in all .m files, otherwise errors will result.
- Press Run and follow instructions in the Command Window
- A Test Matrix featuring the test parameters should accompany the script to enable accurate filenames and sub-folders to be created
- The codes are to be used as a guide only and are to be adjusted depending on the user requirements, e.g. ADV sampling frequency, flume geometry, measurement coordinates, traverse speed, data windowing, filepaths, etc.


NOTES:
- Vertical and cross-stream measurement profiles contain a measurement increment of 100 mm (see 'cruciform_measurement_coordinates.png')
- Vertical profiles contain nine measurement points throughout the 1000 mm channel depth (see 'cruciform_measurement_coordinates.png')
- Horitzontal profiles contain ten measurement points across the 1100 mm channel width (see 'cruciform_measurement_coordinates.png')
- Transect measurement increment is 125 mm in both the horizontal and vertical directions (see 'cruciform_measurement_coordinates.png')
- Measurements performed on a 2D traverse, the speed and acceleration parameters for which can be modified in profiling_post_processing.m
- The folder 'reference' contains an example text matrix with pre-populated inputs (a filename convention for reading the input files is already defined for the script)


TO FOLLOW:
- A script for ADV data analysis and its associated functions will follow at a later date
