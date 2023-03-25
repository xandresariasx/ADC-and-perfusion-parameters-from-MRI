
clear all
Patient='101-003-108';
Options.UseNetworkMasks=0;
    
Options.DoDWI=1;
Options.CurateDWI=0;
    Options.BiasCorrectDWI=0;
Options.AnnotateDWI=1;
Options.GlobalRegisterDWI=0;
Options.LocalRegisterDWI=0;
    Options.roiFactorDW=0;
Options.GetParDWI=0;
Options.GlobalDateRegisterDWI=0;
Options.LocalDateRegisterDWI=0;
    Options.roiFactorLDW=3;
    Options.UseSegmentationForLocalRegDWI=1;
Options.DoT1W=0;
Options.CurateT1W=0;
    Options.BiasCorrectT1=0;
Options.AnnotateT1W=0;
Options.GlobalRegisterT1W=0;
Options.SpetialRegT1=0;
Options.LocalRegisterT1W=0;
    Options.roiFactorT1=[4];                     % 3;
    Options.IncludeArteryLR=[0];
Options.GetT1=0;
    Options.AirTh=0.01;                         % <0.035
    Options.ScaleT1Ims=0;
    Options.RefTiss='Muscle';
    Options.gammaT1 = 10*5e-5;                     % 5e-5;
    Options.gammaM0 = 10*5e-5;                     % 5e-5;
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
    Options.roiFactorLT1=2;                    % 4 
    Options.UseSegmentationForLocalRegT1=1;     % 1
Options.DoT1WDWIReg=0;
    Options.roiFactorLT1DW=2;
    Options.RotateViewDW=[];             % [],'AxToCor'

disp(Patient)
MainP3(Patient, Options)

