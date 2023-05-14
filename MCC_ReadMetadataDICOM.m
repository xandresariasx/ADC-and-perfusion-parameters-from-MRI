% MCC_ReadMetadata  Read metadata from a DICOM file.
%
%   metadata = MCC_ReadMetadataDICOM(filename)
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center


function metadata = MCC_ReadMetadataDICOM(filename)
try
    metadata =  dicominfo(filename);
catch
    CD = cd;
    cd('C:\Program Files\MATLAB\R2015b\toolbox\images\iptformats');
    metadata =  MCC_dicominfo(filename);
    cd(CD);
end
