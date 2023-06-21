# ADC-and-perfusion-parameters-from-MRI
These scripts are intended to reproduce the results of the paper: Arias et al. "Magnetic Resonance Imaging of Tumor Microenvironmental Response to Pegvorhyaluronidase alpha (PEGPH20) in Patients with Advanced Solid Tumors".

# Requirements: 
This tool was developed and tested in Matlab R2017b on Windows OS. It is not guaranteed this tool is going to work on different settings.

# Installation: 
Code: First download all files and add it to the Matlab path. Second, edit RunOnePatient.m. Change the raw data folder, the folder to write the data, and the patient to process. 
Change other options if desired.

Data: The anonymized MRI images used in this study are expected to be shared soon in Imaging Data Commons (https://datacommons.cancer.gov/repository/imaging-data-commons).

# Results: 
For DW-MRI data:

-Curated data is going to be located at: {WriteFolder}\ {Patient}\DWI\ {Date}\

-Registered images at: {WriteFolder}\ {Patient}\DWI\Registered\ {Date}\

-ADC maps at: {WriteFolder}\ {Patient}\DWI\Registered\ {Date}\Local\ADC.mat

-Registered ADC across dates at: {WriteFolder}\ {Patient}\DWI\DateRegisteredLocal\ {Date}\Processed\Registered_VOI#\ADC\

# Copyrights: 
If part of this code is used in your work, please cite:

Arias et al. "Magnetic Resonance Imaging of Tumor Microenvironmental Response to Pegvorhyaluronidase alpha (PEGPH20) in Patients with Advanced Solid Tumors"

If the T1 mapping is also used, also cite:

Sandmann Constantin, Hodneland Erlend,Modersitzki Jan. A practical guideline for T1 reconstruction from various flip angles in MRI Journal of Algorithms & Computational
Technology. 2016
