% MCC_eType Data type enumerator (orginal, resampled, registered ...)
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

classdef MCC_eType
   enumeration
      Original,
      Resampled, 
      Registered,  
      Registered_Elastic,
      Registered_VOI1,
      Registered_VOI2,
      Registered_VOI3,
      Registered_VOI4,
      Registered_VOI5,   % 8/9/2020
      Registered_VOI6,   % 8/9/2020
      Registered_VOI7,   % 8/9/2020
      Registered_VOI8,   % 8/9/2020
      Registered_Date,
      Ready,
      Segmented,
      Temporary,
      Localized      
      Mask,      
      Parametric,
      Parametric_VOI1,
      Parametric_VOI2,
      Quantified,
      DICOM
      All
   end
end