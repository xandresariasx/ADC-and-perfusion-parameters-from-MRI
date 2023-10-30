function VisualizeLocalRegistrationDWI(Folder,LocalRegTemplates,LocalRegTemplate,roiFactor)
% Function to perform local registration around tumor in DW-MRI data

try, rmdir([Folder filesep 'Processed' filesep], 's'); end
try, rmdir([Folder filesep 'Tmp' filesep], 's'); end
ImagesR=AdjustDirVariable(dir(Folder)); 
ImagesR=ImagesR([ImagesR(:).isdir]);
aux=strcmp({ImagesR(:).name},'Processed');
ImagesR(aux)=[];
aux=strsplit(Folder,filesep);
aux(cellfun(@isempty, aux))=[];
Folder2=strjoin(aux(1:end-1),filesep);
try, rmdir([Folder2 filesep 'Local' filesep], 's'); end

Names={ImagesR(:).name};
if roiFactor==0   
    cellfun(@(x) mkdir([Folder2 filesep 'Local' filesep x]),...
        Names)
    cellfun(@(x) copyfile([Folder x], [Folder2 filesep 'Local' filesep x]),...
        Names)
    return
end

aux=cellfun(@(x) strsplit(x,'_'), Names,'UniformOutput', false);
aux2=cellfun(@(x) x(end), aux);
mkdir([Folder 'Tmp' filesep])
for I=unique(aux2)
    I=I{1};
    cellfun(@(x) mkdir([Folder 'Tmp' filesep I filesep x]),...
        Names(strcmp(aux2,I)))
    cellfun(@(x) copyfile([Folder x], [Folder 'Tmp' filesep I filesep x]),...
        Names(strcmp(aux2,I)))
    cellfun(@(x) mkdir([Folder 'Tmp' filesep I filesep 'Processed' filesep 'Registered' filesep x]),...
        Names(strcmp(aux2,I)))
    cellfun(@(x) copyfile([Folder x], [Folder 'Tmp' filesep I filesep 'Processed' filesep 'Registered' filesep x]),...
        Names(strcmp(aux2,I)))
end

for I=unique(aux2)
    I=I{1};
    aux3=Names(strcmp(aux2,I));
    aux4=cellfun(@(x) strsplit(x,'_'), aux3,'UniformOutput', false);
    aux4=cellfun(@(x) [x{1}(3:end) '_' x{2}],aux4,'UniformOutput', false);
    LocalRegTemplates(str2num(I))=aux3(strcmp(aux4,LocalRegTemplates{str2num(I)}));
end

[~,Infos]=ReadDcmFolder4([Folder Names{1} filesep]);

VOI=GenerateVOIforLocalDateRegistration([],...
        Folder,...
        [], [],roiFactor,Infos{1}{1});
    
obj = MCC_cPatient([Folder 'Tmp' filesep]);

obj.Register_Date_Local(LocalRegTemplates,...
    aux2(strcmp(Names,LocalRegTemplate)),...
            {'none','rigid','similarity','affine'},VOI);
        
 for I=unique(aux2)
    I=I{1};
    cellfun(@(x) mkdir([Folder2 filesep 'Local' filesep x]),...
        Names(strcmp(aux2,I)))
    cellfun(@(x) copyfile([Folder 'Tmp' filesep I filesep 'Processed' filesep 'Registered_VOI' num2str(numel(VOI)) filesep x],...
        [Folder2 filesep 'Local' filesep x]),...
        Names(strcmp(aux2,I)))
 end
        
%try, rmdir([Folder '\Tmp\'], 's'); end
 
VolsGR=cellfun(@(x) ReadDcmFolder4([Folder2 filesep 'Global' filesep x filesep]), Names);
[VolsR,InfosLR]=cellfun(@(x) ReadDcmFolder4([Folder2 filesep 'Local' filesep x filesep]), Names);

for I=1:length(VolsR)
    VolsR{I}(VolsR{I}==0)=VolsGR{I}(VolsR{I}==0);
    for J=1:size(VolsR{I},3)
        dicomwrite(uint32(VolsR{I}(:,:,J)),...
            InfosLR{I}{J}.Filename,...
            InfosLR{I}{J}, 'MultiframeSingleFile',false,'CreateMode','copy');
    end
end    
       
disp('Finished local registration')

        
        
        
        

