%MCC_RemoveSeriesNumber  Remove series number label ([ S# ]).
%
%   r = MCC_RemoveSeriesNumber(s)
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

function r = MCC_RemoveSeriesNumber(s)

ind = strfind(s, ' [ S#'); % Image mode case (S# = series number)
if ~isempty(ind)
    s = s(1:(ind-1));
end
r = s;
