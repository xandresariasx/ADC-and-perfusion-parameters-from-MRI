%MCC_GetInstituteSignature Get institute signature.
%
%   r = MCC_GetInstituteSignature(content, pattern)
%
%   INPUTS
%       ~
%
%   OUTPUT
%       r : author signature
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

% Versions:
%
% v134
%   NGR 2016 04 01
%       - Creation

function r = MCC_GetInstituteSignature()

c = clock;

r = [num2str(c(1)) ' © Moffitt Cancer Center'];

end