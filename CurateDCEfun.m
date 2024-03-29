function CurateDCEfun(Infos,SeriesDescription,PatientName,Date,Folder_DCE,WriteFolder)

clc
disp(['Curating DCE Patient: ' PatientName ' Date: ' Date])
CombineAcqusitionAndTriggerTimes=1;

if  any(cellfun(@(x) strcmp(PatientName,x),...
            {'102-003-104'}))
    CombineAcqusitionAndTriggerTimes=0;
end


aux2=cellfun(@(x) ~isempty(regexpi(x,Folder_DCE)),...
        SeriesDescription, 'UniformOutput', false); 
    
aux2=cell2mat(aux2);

InfosDCM=Infos(aux2);

if isempty(InfosDCM)
    return
end

if CombineAcqusitionAndTriggerTimes
    TrigTimes=cell2mat(cellfun(@GetRealTime2, InfosDCM, 'UniformOutput', false));
else
    TrigTimes=cellfun(@(x) x.AcquisitionTime, InfosDCM, 'UniformOutput', false);
    TrigTimes=cell2mat(cellfun(@(x) GetRealTime(x), TrigTimes, 'UniformOutput', false));
end
        
[UTimes,~,Ib]=unique(TrigTimes);

DUTimes=UTimes(2:end)-UTimes(1:end-1);
aux=find(DUTimes<=1);
for I=1:length(aux)
    Ib(Ib==aux(I)+1)=aux(I);
end
InfosPerTime=[];
for I=unique(Ib)'
    aux3=I==Ib;
    InfosPerTime{end+1}=InfosDCM(aux3);
end

%%%%%
aux=cellfun(@length, InfosPerTime);
InfosSeparatedAll=cell(0);
for I=1:length(InfosPerTime)    
    Times=cellfun(@GetTriggerTime , InfosPerTime{I}, 'UniformOutput', false);
    Times=cell2mat(Times);
    [Utimes2,~,Ib2]=unique(Times);
    InfosSeparated=cell(0);
    for J=1:max(Ib2)
        InfosSeparated{J}=InfosPerTime{I}(Ib2==J);        
    end
    InfosSeparatedAll{I}=InfosSeparated;
end

InfosPerTime2=[];
for I=1:length(InfosPerTime)
    InfosPerTime2=[InfosPerTime2 InfosSeparatedAll{I}];
end
InfosPerTime=InfosPerTime2;
%%%%%


for J=1:length(InfosPerTime)
    aux4=[];
    for I=1:length(InfosPerTime{J})
        aux4(I)=InfosPerTime{J}{I}.InstanceNumber;
    end
    if numel(unique(aux4))~=numel(aux4)     % Check
        [~,aux]=unique(aux4);
        aux4=aux4(aux);
    end
    [~,ia]=sort(aux4);
    InfosPerTimeSorted{J}={InfosPerTime{J}{ia}};
    Vols{J}=StackImage(cellfun(@(x) dicomread(x.Filename), InfosPerTimeSorted{J},'UniformOutput' , false));    
end

%%% New
if CombineAcqusitionAndTriggerTimes
    Times=cellfun(@(x) GetRealTime2(x{1}), InfosPerTimeSorted);
else
    Times=cellfun(@(x) x{1}.AcquisitionTime, InfosPerTimeSorted, 'UniformOutput', false);
    Times=cell2mat(cellfun(@(x) GetRealTime(x), Times, 'UniformOutput', false));
end
Times=Times-Times(1);
% if strcmp(PatientName,'102-003-105')
%     any(cellfun(@(x) strcmp(PatientName,x),...
%             {'102-003-105','102-003-117'}))
if Times(10)-Times(9)>1000
    Times=Times/1000;
end
%%%
if exist([WriteFolder PatientName filesep 'T1W' filesep Date filesep])==0
    mkdir([WriteFolder PatientName filesep 'T1W' filesep Date filesep]);
end


aux=combnk([1:numel(Vols)],2);
RM=[];
for I=1:size(aux,1)
    if isequal(Vols{aux(I,1)},Vols{aux(I,2)})
        RM=[RM max(aux(I,1),aux(I,2))];
    end
end
if ~isempty(RM)
    disp(['Automatically discarded volumes: ' num2str(RM)])
end
kdisc=unique([RM]);

fid=fopen([WriteFolder PatientName filesep 'T1W' filesep Date filesep 'DiscardedDCE.txt'],'w');
fprintf(fid, num2str(kdisc));
fclose(fid);


Vols(kdisc)=[];
InfosPerTimeSorted(kdisc)=[];

xdims=cell2mat(cellfun(@(x) size(x,1), Vols,'UniformOutput' , false)');
ydims=cell2mat(cellfun(@(x) size(x,2), Vols,'UniformOutput' , false)');
NSlices=cell2mat(cellfun(@(x) size(x,3), Vols,'UniformOutput' , false)');
aux1=cell2mat(cellfun(@(x) size(x,3)==mode(NSlices), Vols,'UniformOutput' , false));
aux2=cell2mat(cellfun(@(x) size(x,1)==mode(xdims), Vols,'UniformOutput' , false));
aux3=cell2mat(cellfun(@(x) size(x,2)==mode(ydims), Vols,'UniformOutput' , false));
aux=aux1 & aux2 & aux3;
Vols2=Vols(aux);
Vols4D(1,1,1,1:length(Vols2))=Vols2; 
Volav=sum(cell2mat(Vols4D),4)/length(Vols2);

 for k=1:length(InfosPerTimeSorted)
    mkdir([WriteFolder PatientName filesep 'T1W' filesep Date filesep 'DCE_t=' num2str(k) filesep])
    for I=1:length(InfosPerTimeSorted{k})
        copyfile(InfosPerTimeSorted{k}{I}.Filename,[WriteFolder PatientName filesep 'T1W' filesep Date filesep 'DCE_t=' num2str(k) filesep num2str(I) '.dcm'])
    end   
 end


WriteDicomFolderV3(Volav, InfosPerTimeSorted{find(aux,1)}, [WriteFolder PatientName filesep 'T1W' filesep Date filesep 'DCEav' filesep],'DCE Average');

close all force




function Out=ContrastPresent(Info)
Out=0;
try
    aux=Info.ContrastBolusAgent;
    if ~isempty(aux)
        Out=1;
    end
end





