% MCC_WindowButtonDownFcn Slice navigation using graphical input from mouse
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function MCC_WindowButtonDownFcn(src , callbackdata)

global eViewAxis

% Axis
a = gca;
switch a.Title.String
    case char(MCC_eView.Axial)
        eViewAxis = MCC_eView.Axial;
    case char(MCC_eView.Coronal)
        eViewAxis = MCC_eView.Coronal;
    case char(MCC_eView.Sagittal)
        eViewAxis = MCC_eView.Sagittal;
end

end