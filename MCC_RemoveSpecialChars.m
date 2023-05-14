%MCC_RemoveSpecialChars  Remove special charaters to avoid in folder names.
%
%   r = MCC_RemoveSpecialChars(s)
%
%   INPUTS
%       s : string
%
%   OUTPUT
%       r : string
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function r = MCC_RemoveSpecialChars(s)

% Characters to avoid in folder name
ca = {'''','<','>',':','"','?','/','\','|','?','*'};
r = s;
for iCA = 1:length(ca)    
        r = strrep(r, ca(iCA), ' ');    
end