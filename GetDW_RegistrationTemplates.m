function [Templates,TemplateDate,LocalRegTemplate]=GetDW_RegistrationTemplates(Folder)

Foldersb=AdjustDirVariable(dir(Folder));
Foldersb=Foldersb([Foldersb(:).isdir]);
aux=arrayfun(@(x) AdjustDirVariable(dir([x.folder filesep x.name])), Foldersb,'UniformOutput',false);  
if isempty(aux)
    Templates=[];TemplateDate=[];LocalRegTemplate=[];
    return
end
Templates=cellfun(@SelectMinB, aux,'UniformOutput',false);
aux2=sort({Foldersb(:).name});
res=1;
TemplateDate=aux2{res};
LocalRegTemplate=['b=' Templates{res} '_' TemplateDate(3:end)];




function MinB=SelectMinB(Dir)

aux2={Dir.name};
aux3=strcmp(aux2,'Processed');
aux2(aux3)=[];
MinB=sort(aux2);
MinB=MinB{1};

    




