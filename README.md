# ADC-and-perfusion-parameters-from-MRI
These scripts are intended to reproduce the results of the paper: Arias et al. "Magnetic Resonance Imaging of Tumor Microenvironmental Response to Pegvorhyaluronidase alpha (PEGPH20) in Patients with Advanced Solid Tumors". It generates ADC maps from DW-MRI data; T1 maps from T1-MRI; and Ktrans; Vp, Ve, iAUC maps from DCE-MRI.

# Requirements: 
Matlab R2017b or above.

# Installation: 
Code: Download all files and add it to the Matlab path. 

Data: The anonymized MRI images used in this study are expected to be shared soon in Imaging Data Commons (https://datacommons.cancer.gov/repository/imaging-data-commons). Temporary, one patient data with ADC images can be downloaded here: https://drive.google.com/file/d/1N7FeGAPMVt-74yc4petkc5WII1sV-1up/view?usp=drive_link, another patient data with T1 and DCE images here: https://drive.google.com/file/d/1UfuuWHuaUZlcfY0uRNUHsdQgUXO8EDdV/view?usp=sharing. Extract the files in the raw data folder.

# Usage:
Edit RunOnePatient.m. Change the raw data folder, the folder to write the data, and the patient to process. 
Change other options if desired. Then execute RunOnePatient.m.


# Results: 
For DW-MRI data:

-Curated data is going to be located at: {WriteFolder}\{Patient}\DWI\{Date}\

-Registered images at: {WriteFolder}\{Patient}\DWI\Registered\{Date}\

-ADC maps at: {WriteFolder}\{Patient}\DWI\Registered\{Date}\Local\ADC.mat

-Registered ADC across dates at: {WriteFolder}\{Patient}\DWI\DateRegisteredLocal\{Date}\Processed\Registered_VOI#\ADC\

For T1-MRI and DCE-MRI data:

-Curated data is going to be located at: {WriteFolder}\{Patient}\T1W\{Date}\

-Registered images at: {WriteFolder}\{Patient}\T1W\Registered\{Date}\Local\

- T1 maps at: {WriteFolder}\{Patient}\T1W\Registered\{Date}\Local\T1Data_NoCorrection.mat

- iAUC maps at: {WriteFolder}\{Patient}\T1W\Registered\{Date}\Local\90sAUC.mat

- Ktrans, Vp, Ve maps at: {WriteFolder}\{Patient}\T1W\Registered\{Date}\Local\Perfussion_Parameters_Maps.mat



# Copyrights: 
If part of this code is used in your work, please cite:

Arias et al. "Magnetic Resonance Imaging of Tumor Microenvironmental Response to Pegvorhyaluronidase alpha (PEGPH20) in Patients with Advanced Solid Tumors"

If the T1 mapping is also used, also cite:

Sandmann Constantin, Hodneland Erlend,Modersitzki Jan. A practical guideline for T1 reconstruction from various flip angles in MRI Journal of Algorithms & Computational
Technology. 2016
