function OrganizeRegisteredDWIP3(Folder,WriteFolder)


ImagesDate=AdjustDirVariable(dir(Folder));

aux=strsplit(Folder,'\'); 
FolderD=[strjoin(aux(1:end-1),'\') '\Registered\' aux{end} '\Global\'];

for Image1={ImagesDate.name}
    Images=AdjustDirVariable(dir([Folder '\' Image1{1} '\Processed\Registered_Date\']));    
    for Image={Images.name}        
        aux3=strsplit(Image1{1},'='); 
        aux2=strsplit(Image{1},'_');
        aux4=aux2{1};
        copyfile([Folder '\' Image1{1} '\Processed\Registered_Date\' Image{1} '\'],...
                [FolderD 'b=' Image{1} '_' aux3{end} '\'])         
    end
end