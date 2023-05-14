%MCC_Disp  Display value of variable with specified tabulation
%
%   ~ = MCC_Disp(txt, Ntab)
%
%   INPUTS
%        txt: string to display
%       Ntab: number of tabulations
%
%   OUTPUT
%       ~
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function MCC_Disp(txt, Ntab)

if nargin<=1
    Ntab = 2;
end
for i=1:Ntab
    txt = [sprintf('\t') txt];
end
disp(txt);