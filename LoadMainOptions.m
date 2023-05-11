
function [DoDWI,CurateDWI,BiasCorrectDWI,AnnotateDWI,...
GlobalRegisterDWI,LocalRegisterDWI,roiFactorDW,GetParDWI,GlobalDateRegisterDWI,...
LocalDateRegisterDWI,roiFactorLDW,UseSegmentationForLocalRegDWI,DoT1W,...
CurateT1W,BiasCorrectT1,AnnotateT1W,GlobalRegisterT1W,SpetialRegT1,LocalRegisterT1W,...
roiFactorT1,IncludeArteryLR,GetT1,CorrectFA,WeightT1,AirTh,ScaleT1Ims,RefTiss,gammaT1,gammaM0,...
solver,CorrectBiasBetweenSlices,AssumeT1s,GetConc,ContrastBrand,ScaleMo,AssumeT1sC,AssumeT1Artery,GetPar,...
GlobalDateRegisterT1,LocalDateRegisterT1,roiFactorLT1,UseSegmentationForLocalRegT1,...
DoT1WDWIReg,roiFactorLT1DW,RotateViewDW]=LoadMainOptions(Options)
% Function that loads user options    
    
DoDWI=Options.DoDWI;
CurateDWI=Options.CurateDWI;
    BiasCorrectDWI=Options.BiasCorrectDWI;
AnnotateDWI=Options.AnnotateDWI;
GlobalRegisterDWI=Options.GlobalRegisterDWI;
LocalRegisterDWI=Options.LocalRegisterDWI;
    roiFactorDW=Options.roiFactorDW;
GetParDWI=Options.GetParDWI;
GlobalDateRegisterDWI=Options.GlobalDateRegisterDWI;
LocalDateRegisterDWI=Options.LocalDateRegisterDWI;
    roiFactorLDW=Options.roiFactorLDW;
    UseSegmentationForLocalRegDWI=Options.UseSegmentationForLocalRegDWI;
DoT1W=Options.DoT1W;
CurateT1W=Options.CurateT1W;
    BiasCorrectT1=Options.BiasCorrectT1;
AnnotateT1W=Options.AnnotateT1W;
GlobalRegisterT1W=Options.GlobalRegisterT1W;
SpetialRegT1=Options.SpetialRegT1;
LocalRegisterT1W=Options.LocalRegisterT1W;
    roiFactorT1=Options.roiFactorT1;                     
    IncludeArteryLR=Options.IncludeArteryLR;
GetT1=Options.GetT1;
    AirTh=Options.AirTh;                         
    ScaleT1Ims=Options.ScaleT1Ims;
    RefTiss=Options.RefTiss;
    gammaT1 = Options.gammaT1;                     
    gammaM0 = Options.gammaM0;                     
    solver    = Options.solver;           
    CorrectBiasBetweenSlices=Options.CorrectBiasBetweenSlices;
    AssumeT1s=Options.AssumeT1s;
    CorrectFA=Options.CorrectFA;
    WeightT1=Options.WeightT1;
GetConc=Options.GetConc;
    ContrastBrand=Options.ContrastBrand;
    ScaleMo=Options.ScaleMo;
    AssumeT1sC=Options.AssumeT1sC;
    AssumeT1Artery=Options.AssumeT1Artery;    
GetPar=Options.GetPar;
GlobalDateRegisterT1=Options.GlobalDateRegisterT1;
LocalDateRegisterT1=Options.LocalDateRegisterT1;
    roiFactorLT1=Options.roiFactorLT1;                      
    UseSegmentationForLocalRegT1=Options.UseSegmentationForLocalRegT1;     
DoT1WDWIReg=Options.DoT1WDWIReg;
    roiFactorLT1DW=Options.roiFactorLT1DW;
    RotateViewDW=Options.RotateViewDW;    
