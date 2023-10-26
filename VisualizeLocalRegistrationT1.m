function SlicePos=VisualizeLocalRegistrationT1(FolderAux, ImageTemplate, roiFactor, IncludeArtery)


ImagesR=AdjustDirVariable(dir(FolderAux)); 

ImagesNames={ImagesR(:).name};
Ind0=[ImagesR(:).isdir] & ~cell2mat(cellfun(@(x) strcmp(x,'Copy'), ImagesNames,'UniformOutput' , false));
ImagesNames=ImagesNames(Ind0);
ImagesR=ImagesR(Ind0);
[~,Ind]=sort(cell2mat(cellfun(@(x) str2num(x(strfind(x,'t=')+2:end)), ImagesNames,'UniformOutput' , false)));
ImagesRDCE=ImagesR(Ind);

ImagesRT1=ImagesR(cell2mat(cellfun(@(x) contains(x,'a='), ImagesNames,'UniformOutput' , false)));

VolsRDCE=arrayfun(@(x) ReadDcmFolder4([x.folder filesep x.name filesep]), ImagesRDCE);

VolsRT1=arrayfun(@(x) ReadDcmFolder4([x.folder filesep x.name filesep]), ImagesRT1);

ImagesAux={ImagesRDCE(:).name ImagesRT1(:).name};


aux=strsplit(FolderAux,filesep);
aux(cellfun(@isempty, aux))=[];
FolderT1=strjoin(aux(1:find(strcmp(aux,'T1W'))),filesep);
Date=aux{find(strcmp(aux,'T1W'))+1};

DirM=dir([FolderT1 filesep Date filesep 'Processed' filesep 'Registered_VOI*']);
try, cellfun(@(x,y) rmdir([x filesep y filesep],'s'),{DirM(:).folder},{DirM(:).name});end
try, rmdir([FolderT1 filesep 'Registered' filesep Date filesep 'Global' filesep],'s');  end
try, rmdir([FolderT1 filesep 'Registered' filesep Date filesep 'Local' filesep],'s');  end
mkdir([FolderT1 filesep 'Registered' filesep Date filesep 'Global'])
mkdir([FolderT1 filesep 'Registered' filesep Date filesep 'Local'])
cellfun(@(x) copyfile([FolderAux filesep x],[FolderT1 filesep 'Registered' filesep Date filesep 'Global' filesep x]), ImagesAux,'UniformOutput' , false);
cellfun(@(x) copyfile([FolderAux filesep x],[FolderT1 filesep 'Registered' filesep Date filesep 'Local' filesep x]), ImagesAux,'UniformOutput' , false);

if  all(roiFactor==0) && IncludeArtery(1)==0        % 1/28/2021 if no local registration required finish here
    return
end

kdisc=[1:length(ImagesAux)];
ImagesLR=ImagesNames;

obj = MCC_cPatient(FolderT1);
aux2={obj.content{:}};
aux2=~strcmp(aux2,Date);
obj.content(aux2)=[];
obj.StudyDate(aux2)=[];
IndLR=cell2mat(cellfun(@(x) find(strcmp(x,obj.StudyDate.content)), ImagesLR,'UniformOutput',false));
aux=ones(size(obj.StudyDate.content));
aux(IndLR)=0;
obj.StudyDate.content(logical(aux))=[];   
%
aux=AdjustDirVariable(dir([FolderT1 filesep Date filesep 'Processed' filesep]));
aux=aux([aux(:).isdir]);
aux=aux(contains({aux(:).name}, 'VOI'));
for I=1:length(aux)
    rmdir([aux(I).folder filesep aux(I).name filesep], 's')
end
%
load([FolderT1 filesep Date filesep 'Masks.mat'])
MatLabels=cat(4,MaksPerLabel{:});
aux2=[];
if ~isempty(VolsRDCE)
    for I=1:length(MaksPerLabel),aux2{I}=I*ones(size(VolsRDCE{1})); end
else
    for I=1:length(MaksPerLabel),aux2{I}=I*ones(size(VolsRT1{1})); end
end
aux2=cat(4,aux2{:});
MatLabels=MatLabels.*aux2;

obj.StudyDate.LoadData;
obj=GenerateROIsLocalReg(obj,obj.StudyDate.vsRegistered(1).metadata{1}, Labels, MatLabels, IncludeArtery, roiFactor);

for I=1:length(obj.StudyDate(1).VOI)
    try
    obj.Register({ImageTemplate},[],[],['voi' num2str(I)]);
    catch ME
        ME
        disp('local reg. error')
        I=I-1;
        break;
    end
end

ImagesC=obj.StudyDate.content;
cellfun(@(x) copyfile([FolderT1 filesep Date filesep 'Processed' filesep 'Registered_VOI' num2str(I) filesep x],...
    [FolderT1 filesep 'Registered' filesep Date filesep 'Local' filesep x]), ImagesC,'UniformOutput' , false);

%%% NEW - FILL REGISTRATION GAPS, BETTER VISUALLY APPEALLING RESULTS
VolsGR=cellfun(@(x) ReadDcmFolder4([FolderT1 filesep 'Registered' filesep Date filesep 'Global' filesep x filesep]), ImagesC);
[VolsLR,InfosLR]=cellfun(@(x) ReadDcmFolder4([FolderT1 filesep 'Registered' filesep Date filesep 'Local' filesep x filesep]), ImagesC);

for I=1:length(VolsLR)
    VolsLR{I}(VolsLR{I}==0)=VolsGR{I}(VolsLR{I}==0);
    for J=1:size(VolsLR{I},3)
        dicomwrite(uint32(VolsLR{I}(:,:,J)),...
            InfosLR{I}{J}.Filename,...
            InfosLR{I}{J}, 'MultiframeSingleFile',false,'CreateMode','copy');
    end
end    










