%MCC_GetVersionFormated() Get formated version number (#.#.#)
%
%   INPUTS
%       ~
%
%   OUTPUT
%       r : string formated (#.#.#)
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function r = MCC_GetVersionFormated()

v = num2str(MCC_GetVersion);

r = [v(1) '.' v(2) '.' v(3)];

end