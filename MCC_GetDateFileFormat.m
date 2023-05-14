%MCC_GetDateFileFormat  Get data file format.
%
%   r = MCC_GetDateFileFormat()
%
%   INPUTS
%       ~
%
%   OUTPUT
%       r : formated date (yyyy_mm_dd)
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function r = MCC_GetDateFileFormat()

c = clock;

y = num2str(c(1));

if c(2)<=9
    m = ['0' num2str(c(2))];
else
    m = num2str(c(2));
end

if c(3)<=9
    d = ['0' num2str(c(3))];
else
    d = num2str(c(3));
end

r = [y '_' m '_' d];

end