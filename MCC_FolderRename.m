%MCC_FolderRename  Rename all folders.
%
%   MCC_FolderRename(root)
%
%   Radiomics project: "Multi-parametric MRI (mpMRI) Analysis"%
%   Nicolas Georges Rognin, PhD
%   2016-2017 © Moffitt Cancer Center

function MCC_FolderRename(root)

nameFolds = MCC_FolderList(root);

CD = cd;
cd(root);

for i=1:length(nameFolds)
    str = nameFolds{i};
    idx  = strfind(str, ' [ S#');
    if ~isempty(idx)    
        str_new = str(1:(idx-1));
        cmd = ['!rename "' str ' " "' str_new '"'];        
        eval(cmd);
    end
end

cd(CD);