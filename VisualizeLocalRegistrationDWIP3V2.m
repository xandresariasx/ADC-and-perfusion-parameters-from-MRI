function VisualizeLocalRegistrationDWIP3V2(Folder,LocalRegTemplates,LocalRegTemplate,roiFactor)


try, rmdir([Folder '\Processed\'], 's'); end
try, rmdir([Folder '\Tmp\'], 's'); end
ImagesR=AdjustDirVariable(dir(Folder)); 
ImagesR=ImagesR([ImagesR(:).isdir]);
aux=strcmp({ImagesR(:).name},'Processed');
ImagesR(aux)=[];
aux=strsplit(Folder,'\');
aux(cellfun(@isempty, aux))=[];
Folder2=strjoin(aux(1:end-1),'\');
try, rmdir([Folder2 '\Local\'], 's'); end

Names={ImagesR(:).name};
if roiFactor==0   % 1/21/2022 in case roi=0
    cellfun(@(x) copyfile([Folder x], [Folder2 '\Local\' x]),...
        Names)
    return
end

aux=cellfun(@(x) strsplit(x,'_'), Names,'UniformOutput', false);
aux2=cellfun(@(x) x(end), aux);
mkdir([Folder 'Tmp\'])
for I=unique(aux2)
    I=I{1};
    cellfun(@(x) copyfile([Folder x], [Folder 'Tmp\' I '\' x]),...
        Names(strcmp(aux2,I)))
    cellfun(@(x) copyfile([Folder x], [Folder 'Tmp\' I '\Processed\Registered\' x]),...
        Names(strcmp(aux2,I)))
end

for I=unique(aux2)
    I=I{1};
    aux3=Names(strcmp(aux2,I));
    aux4=cellfun(@(x) strsplit(x,'_'), aux3,'UniformOutput', false);
    aux4=cellfun(@(x) [x{1}(3:end) '_' x{2}],aux4,'UniformOutput', false);
    LocalRegTemplates(str2num(I))=aux3(strcmp(aux4,LocalRegTemplates{str2num(I)}));
end

[~,Infos]=ReadDcmFolder4([Folder Names{1} '\']);

VOI=GenerateVOIforLocalDateRegistration([],...
        Folder,...
        [], [],roiFactor,Infos{1}{1});
    
obj = MCC_cPatient([Folder 'Tmp\']);

obj.Register_Date_Local(LocalRegTemplates,...
    aux2(strcmp(Names,LocalRegTemplate)),...
            {'none','rigid','similarity','affine'},VOI);
        
 for I=unique(aux2)
    I=I{1};
    cellfun(@(x) copyfile([Folder 'Tmp\' I '\Processed\Registered_VOI' num2str(numel(VOI)) '\' x],...
        [Folder2 '\Local\' x]),...
        Names(strcmp(aux2,I)))
 end
        
%try, rmdir([Folder '\Tmp\'], 's'); end
 
VolsGR=cellfun(@(x) ReadDcmFolder4([Folder2 '\Global\' x '\']), Names);
[VolsR,InfosLR]=cellfun(@(x) ReadDcmFolder4([Folder2 '\Local\' x '\']), Names);

for I=1:length(VolsR)
    VolsR{I}(VolsR{I}==0)=VolsGR{I}(VolsR{I}==0);
    for J=1:size(VolsR{I},3)
        dicomwrite(uint32(VolsR{I}(:,:,J)),...
            InfosLR{I}{J}.Filename,...
            InfosLR{I}{J}, 'MultiframeSingleFile',false,'CreateMode','copy');
    end
end    
       
disp('Finished local registration')

        
        
        
        

