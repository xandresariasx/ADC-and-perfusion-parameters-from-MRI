%MCC_FolderList  Get folder list.
%
%   nameFolds = MCC_FolderList(root)
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"%
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function nameFolds = MCC_FolderList(root)

d = dir(root);
isub = [d(:).isdir]; % returns logical vector
nameFolds = {d(isub).name};
N = length(nameFolds);
if (N<=2)
    nameFolds = [];
else
    nameFolds = nameFolds(3:end);
end