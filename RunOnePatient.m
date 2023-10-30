% Input patient ID and select options to process images
clear all
Folder='Raw/';            % Folder with raw images
WriteFolder='Images/';          % Folder to write processed data
Patient='101-002-114';            % Select patient ID (see paper for available IDs).                                                           
    
Options.DoDWI=1;                             % Select this option to process DW-MRI images.
Options.CurateDWI=0;                       % Option to separate DW images from raw data folder.     
Options.GlobalRegisterDWI=0;           % Option to register the images per visit.
Options.AnnotateDWI=0;                   % Option to manually annotate tissues in DW-MRI images. Name labels as: Tumor1, Tumor2,...                                            
Options.LocalRegisterDWI=0;             % Option to do local registration around tumors. If local registration is not needed the options below should be 0.  
    Options.roiFactorDW=1.5;              % Factor to resize tumor roi (0: no local registration).
Options.GetParDWI=0;                       % Option to get ADC maps                  
Options.GlobalDateRegisterDWI=1;    % Register ADC maps across dates globally.
Options.LocalDateRegisterDWI=1;      % Register ADC maps across dates locally.
    Options.roiFactorLDW=0;                % ROI factor: ROI is equal to factor multiplied by bounding box around tumors.
    Options.UseSegmentationForLocalRegDWI=0;    % If 1, use segmentations instead of image intensities for registration. 
Options.DoT1W=0;                % Option to obtain T1 and DCE parameter maps.
Options.CurateT1W=1;          % Option to separate T1 and DCE images from raw data folder.
Options.AnnotateT1W=0;      % Option to manually annotate tissues in T1-MRI images.
Options.GlobalRegisterT1W=0;    % Option to register the images per visit.
Options.SpetialRegT1=0;             % Option to do registration for special cases that would not work with regular registation such as images with very different FOV.        
Options.LocalRegisterT1W=0;      % Option to do local registration around tumors and artery. If local registration is not needed the options below should be 0.  
    Options.roiFactorT1=[1.5];       % ROI factor per tumor: ROI size is equal to factor multiplied by bounding box around tumors.            
    Options.IncludeArteryLR=[1];   % 1 if the artery should be included in the local registration
Options.GetT1=0;                     % Option to compute T1 maps
    Options.AirTh=0.01;             % Ratio threshold for air in image ( for instance 0.01 means any voxel below 10% of max intensity is air)                      
    Options.ScaleT1Ims=0;         % If different image scaling between T1 images (requires curve fitting toolbox).
        Options.RefTiss='Muscle';  % Input reference tissue in case ScaleT1Ims=1;
    Options.gammaT1 = 10*5e-5;    % Regularization weight for T1 maps.                   
    Options.gammaM0 = 10*5e-5;    % Regularization weight for M0 maps.                   
    Options.solver    = 'pcg';            % 'MatlabInternal'  'pcg';
    Options.WeightT1=1;         % 1 (no weighting to images to get T1 maps), 2 (weight is inverse to noise), 3 (weight is inverse to noise and proportional to max intensity)
    Options.AssumeT1s=0;       % Assume T1 maps.
    Options.CorrectFA=1;         % Correct FA per image in case of big intensity bias (requires curve fitting toolbox).
Options.GetConc=0;               % Option to compute contrast curves
    Options.ContrastBrand='Magnevist';        % Multihance, Magnevist
    Options.ScaleMo=1;                                 % If different image scaling between DCE and T1-MRI.
    Options.AssumeT1sC=0;                          % Assume T1 maps if real T1 maps not available.
    Options.AssumeT1Artery=0;                    % Assume T1 map at artery only in case T1 at artery is very off due to artifacts.
Options.GetPar=0;                 % Option to compute iAUC, Ktrans, Vp, Ve parameter maps.
Options.GlobalDateRegisterT1=0;         % Register maps across dates globally.
Options.LocalDateRegisterT1=0;           % Register maps acress dates locally.
    Options.roiFactorLT1=2;                    % ROI factor: ROI is equal to factor multiplied by bounding box around tumors.            
    Options.UseSegmentationForLocalRegT1=0;   % If 1, use segmentations instead of image intensities for registration.  
Options.DoT1WDWIReg=0;
    Options.roiFactorLT1DW=2;
    Options.RotateViewDW=[];             

disp(Patient)
Main(Patient, Folder, WriteFolder, Options)

