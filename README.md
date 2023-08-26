# adv_data_processing

There is currently one high level script in the repository for ADV data post-processing (adv_post_process.m).

OVERVIEW:
Scripts reads ASCII file(s) from a Nortek Vectrino and calls multiple functions to post-process the ADV data.

OBJECTIVES:
- Arrange the data into appropriate arrays
- Apply data data windowing
- Apply despiking algorithm (Goring and Nikora (2002))
- Apply rotational corrections to the data to account for any probe misalignment
- Determine noise floor
- Save data

INSTRUCTIONS:
The user is to add their own filepaths in all .m files, otherwise errors will result
Press Run and follow instructions in the Command Window

The linked functions include:
merging.m - 


NOTES:
A script for ADV data analysis and its associated functions will follow at a later date.
