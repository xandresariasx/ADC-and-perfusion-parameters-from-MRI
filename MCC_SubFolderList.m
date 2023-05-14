% MCC_SubFolderList  Returns a the list of sub-folder into a cell array
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function r = MCC_SubFolderList(folder)

lst = dir(folder);
lst(~[lst.isdir])= [];
r = [];
j = 1;
for i=1:length(lst);
    txt = lst(i).name;
    if ~strcmp(txt,'.') && ~strcmp(txt,'..') && ~strcmp(txt,'Processed')
        r{j} = txt;
        j = j +1;
    end
end