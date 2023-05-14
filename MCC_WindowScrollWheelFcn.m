% MCC_WindowScrollWheelFcn Slice navigation using mouse scroll wheel
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function MCC_WindowScrollWheelFcn(src , callbackdata)

try
    
    global eViewAxis
    obj = MCC_getGlobalx;
    
    % Plane
    switch eViewAxis
        case MCC_eView.Axial
            p = 3;
            iP = obj.vPlaneCur(p) + callbackdata.VerticalScrollCount;
        case MCC_eView.Coronal
            p = 1;
            iP = obj.vPlaneCur(p) + callbackdata.VerticalScrollCount;
        case MCC_eView.Sagittal
            p = 2;
            iP = obj.vPlaneCur(p) - callbackdata.VerticalScrollCount;
    end
    V = obj.vsTem(1).data;
    
    if iP>size(V,p)
        iP = size(V,p);
    end
    if iP<1
        iP = 1;
    end
    obj.vPlaneCur(p) = iP;
    vShiftPlane(1) = obj.vPlaneCur(1) - round(size(V,1)/2);
    vShiftPlane(2) = obj.vPlaneCur(2) - round(size(V,2)/2);
    vShiftPlane(3) = obj.vPlaneCur(3) - round(size(V,3)/2);
    vShiftPlane(p) = iP - round(size(V,p)/2);
    
    % Display
    obj.Display4Q(MCC_eType.Original, vShiftPlane, 0, 0);
    
    % Save
    obj.SaveWorkspace;
    
end

end