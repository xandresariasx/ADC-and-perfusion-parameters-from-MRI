% MCC_WindowButtonUpFcn Slice navigation using graphical input from mouse
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center


function MCC_WindowButtonUpFcn(src , callbackdata)

global eViewAxis
obj = MCC_getGlobalx;

% Axis
a = gca;

iVOI = 0;
for i=1:length(obj.hROI)         
    xyzDxDyDz_mm_prev = obj.VOI(i).xyzDxDyDz_mm;    
    switch a.Title.String
        case char(MCC_eView.Axial)
            eViewAxis = MCC_eView.Axial;
            p = obj.hROI(i).h(1).getPosition;
            idx = [1 2 4 5];
            obj.VOI(i).xyzDxDyDz_mm(idx) = p;
        case char(MCC_eView.Coronal)
            eViewAxis = MCC_eView.Coronal;
            p = obj.hROI(i).h(2).getPosition;
            idx = [1 3 4 6];
            obj.VOI(i).xyzDxDyDz_mm(idx) = p;
        case char(MCC_eView.Sagittal)
            p = obj.hROI(i).h(3).getPosition;
            idx = [2 3 5 6];
            obj.VOI(i).xyzDxDyDz_mm(idx) = p;
            eViewAxis = MCC_eView.Sagittal;
    end    
    % Identify selected ROI
    if sum(xyzDxDyDz_mm_prev - obj.VOI(i).xyzDxDyDz_mm)
        iVOI = i;
    end
    obj.VOI_mm(i,:) = obj.VOI(i).xyzDxDyDz_mm;
end

% Data
r = obj.vsTemporary(1).res_mm;
if iVOI
    % Auto-plane selection
    vShiftPlane(2) = round((obj.VOI_mm(iVOI,1) + obj.VOI_mm(iVOI,4)/2)/r(1));
    vShiftPlane(1) = round((obj.VOI_mm(iVOI,2) + obj.VOI_mm(iVOI,5)/2)/r(2));
    vShiftPlane(3) = round((obj.VOI_mm(iVOI,3) + obj.VOI_mm(iVOI,6)/2)/r(3));    
else    
    s = size(obj.vsTemporary(1).data);
    vShiftPlane(1) = obj.vPlaneCur(1) - round(s(1)/2);
    vShiftPlane(2) = obj.vPlaneCur(2) - round(s(2)/2);
    vShiftPlane(3) = obj.vPlaneCur(3) - round(s(3)/2);    
end
obj.Display4Q(obj.eTypeCur, vShiftPlane, 0, 0)

% 3D view
obj.viewpoint = view;

% Save
obj.SaveWorkspace;

end