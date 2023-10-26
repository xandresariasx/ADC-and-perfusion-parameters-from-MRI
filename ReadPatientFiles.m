function [FilesFolderNames,Infos,SeriesDescription]=ReadPatientFiles(FolderBase)
% Read patient files    

Folder=FolderBase;
Files=AdjustDirVariable(dir(Folder));

aux=[Files(:).isdir];
Files(aux)=[];

FileNames={Files(:).name};
FileFolders={Files(:).folder};

FilesFolderNames=cellfun(@(x,y) [x filesep y], FileFolders, FileNames, 'UniformOutput', false);

Isdicoms=cellfun(@(x) isdicom(x), FilesFolderNames, 'UniformOutput', false);

FilesFolderNames(~cell2mat(Isdicoms))=[];;

Infos=cellfun(@dicominfo, FilesFolderNames, 'UniformOutput', false);

try
    SeriesDescription=cellfun(@(x) x.SeriesDescription,Infos, 'UniformOutput', false);
catch
    for I=1:numel(Infos)
        try
            SeriesDescription{I}=Infos{I}.SeriesDescription;
        catch
            SeriesDescription{I}='';
        end    
    end
end



function b=isdicom(x)
try
    dicominfo(x);
    b=1;
catch
    b=0;
end



