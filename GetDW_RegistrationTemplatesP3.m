function [Templates,TemplateDate,LocalRegTemplate]=GetDW_RegistrationTemplatesP3(Folder)

Foldersb=AdjustDirVariable(dir(Folder));
Foldersb=Foldersb([Foldersb(:).isdir]);

aux=arrayfun(@(x) AdjustDirVariable(dir([x.folder '\' x.name])), Foldersb,'UniformOutput',false);  

if isempty(aux)
    Templates=[];TemplateDate=[];LocalRegTemplate=[];
    return
end

Templates=cellfun(@SelectMinB, aux,'UniformOutput',false);

aux2=sort({Foldersb(:).name});

%%%New
close all force
uiopen([Folder '\DWImages.fig'],1)
res=1000;
while res>length({Foldersb(:).name})
    clc
%     aux3=arrayfun(@(x) numel(AdjustDirVariable(dir([x.folder '\'...
%     x.name])))-1, Foldersb);
%     Str=[];
%     for I=1:numel(aux3)
%         Str{I}=['(' num2str(1+(I-1)*aux3(I)) ')' '-' '(' num2str(I*aux3(I)) ')'];
%     end
%     disp(cellfun(@(x,y) [x ' ' y], {Foldersb(:).name}, Str,'UniformOutput', false))
%     res=str2num(input(['Select registration template (default (1): ' aux2{1} ') ? '],'s'));
    Strs=[];
    Str2='Select registration template ';
    for I=1:numel(aux2)
        Str2=[Str2 num2str(I) '. ' aux2{I} ' ('];
        %disp(aux2{I})
        Strs{I}={aux{I}(:).name};
        Strs{I}(cellfun(@(x) strcmp(x,'Processed'), Strs{I}))=[];
        %disp(Strs{I})
        aux4=join(Strs{I});
        Str2=[Str2 aux4{1} '), '];
    end
    Str2=Str2(1:end-2);
    res=str2num(input(Str2,'s'));
end
if isempty(res)
    res=1;
end
%%%

TemplateDate=aux2{res};

LocalRegTemplate=['b=' Templates{res} '_' TemplateDate(3:end)];




function MinB=SelectMinB(Dir)

aux2={Dir.name};
aux3=strcmp(aux2,'Processed');
aux2(aux3)=[];
MinB=sort(aux2);
MinB=MinB{1};

    




