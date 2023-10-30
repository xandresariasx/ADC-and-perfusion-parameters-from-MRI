function OrganizeRegisteredDWI(Folder,WriteFolder)

ImagesDate=AdjustDirVariable(dir(Folder));

aux=strsplit(Folder,filesep); 
FolderD=[strjoin(aux(1:end-1),filesep) filesep 'Registered' filesep aux{end} filesep 'Global' filesep];

for Image1={ImagesDate.name}
    Images=AdjustDirVariable(dir([Folder filesep Image1{1} filesep 'Processed' filesep 'Registered_Date' filesep]));    
    for Image={Images.name}        
        aux3=strsplit(Image1{1},'='); 
        aux2=strsplit(Image{1},'_');
        aux4=aux2{1};
        mkdir([FolderD 'b=' Image{1} '_' aux3{end} filesep])
        copyfile([Folder filesep Image1{1} filesep 'Processed' filesep 'Registered_Date' filesep Image{1} filesep],...
                [FolderD 'b=' Image{1} '_' aux3{end} filesep])         
    end
end
