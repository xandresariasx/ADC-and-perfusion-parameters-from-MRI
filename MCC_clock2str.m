%MCC_clock2str Current date and time as string
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function r = MCC_clock2str()
c = clock;
r = '';
for iC = 1:length(c)
    r = [r ' ' num2str(c(iC))];
end