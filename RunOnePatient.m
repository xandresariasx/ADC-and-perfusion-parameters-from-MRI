% Input patient ID and select options to process images
clear all
Folder='RawImages\';            % Folder with raw images
WriteFolder='Images\';          % Folder to write processed data
Patient='101-002-106';                                                 % Select patient ID (see paper for available IDs)                                                            
    
Options.DoDWI=1;                                                        % Select this option to process DW-MRI images
Options.CurateDWI=1;                                                  %
    Options.BiasCorrectDWI=0;
Options.AnnotateDWI=1;
Options.GlobalRegisterDWI=1;
Options.LocalRegisterDWI=1;
    Options.roiFactorDW=0;
Options.GetParDWI=1;
Options.GlobalDateRegisterDWI=1;
Options.LocalDateRegisterDWI=1;
    Options.roiFactorLDW=3;
    Options.UseSegmentationForLocalRegDWI=1;
Options.DoT1W=1;
Options.CurateT1W=1;
    Options.BiasCorrectT1=0;
Options.AnnotateT1W=1;
Options.GlobalRegisterT1W=1;
Options.SpetialRegT1=0;
Options.LocalRegisterT1W=0;
    Options.roiFactorT1=[4];                   
    Options.IncludeArteryLR=[0];
Options.GetT1=0;
    Options.AirTh=0.01;                         
    Options.ScaleT1Ims=0;
    Options.RefTiss='Muscle';
    Options.gammaT1 = 10*5e-5;                     
    Options.gammaM0 = 10*5e-5;                     
    Options.solver    = 'pcg';            % 'MatlabInternal'  'pcg';
    Options.CorrectBiasBetweenSlices=0;
    Options.AssumeT1s=0;
    Options.CorrectFA=0;
    Options.WeightT1=1;         % 1, 2, 3
Options.GetConc=0;
    Options.ContrastBrand='Magnevist';        % Multihance, Magnevist
    Options.ScaleMo=0;
    Options.AssumeT1sC=0;
    Options.AssumeT1Artery=0;
Options.GetPar=0;
Options.GlobalDateRegisterT1=0;
Options.LocalDateRegisterT1=1;
    Options.roiFactorLT1=2;                    
    Options.UseSegmentationForLocalRegT1=1;     
Options.DoT1WDWIReg=0;
    Options.roiFactorLT1DW=2;
    Options.RotateViewDW=[];             

disp(Patient)
Main(Patient, Folder, WriteFolder, Options)

