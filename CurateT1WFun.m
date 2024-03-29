function CurateT1WFun(WriteFolder,Infos,SeriesDescription,PatientName,Date,Folder_T1) 

clc
close all force
disp(['Curating T1W Patient: ' PatientName ' Date: ' Date])
try
rmdir([WriteFolder PatientName filesep 'T1W' filesep Date],'s')
end
 
aux2=cellfun(@(x) ~isempty(regexpi(x,Folder_T1)),...
        SeriesDescription, 'UniformOutput', false); 

aux2=cell2mat(aux2);

InfosT1s=Infos(aux2);
if isempty(InfosT1s)
    return
end

SeriesDescriptionT1=SeriesDescription(aux2);
[Ut1,~,It1]=unique(SeriesDescriptionT1);

for I=1:length(Ut1)
    InfosT1sI=InfosT1s(I==It1);
    alpha=unique(cellfun(@(x) x.FlipAngle, InfosT1sI));
    try
%         aux=cell2mat(cellfun(@(x) ContrastPresent(x), InfosT1sI, 'UniformOutput', false));
%         InfosT1sI=InfosT1sI(~aux);
    end
    Series=cellfun(@(x) x.SeriesInstanceUID, InfosT1sI, 'UniformOutput', false);
    [SeriesU,~,Ib]=unique(Series);
    if length(SeriesU)<10      % DCE so not use  
        for J=1:length(SeriesU)
            InfosT1sIJ=InfosT1sI(J==Ib);
            if length(InfosT1sIJ)<100
                % Several images with the same serie ID
                SL=cell2mat(cellfun(@(x) x.SliceLocation, InfosT1sIJ, 'UniformOutput', false));
                [a,b,c]=unique(SL);
                for II=1:numel(SL)/max(b)
                    InfosT1sIJ_II{II}=InfosT1sIJ(1+(II-1)*max(b):II*max(b)); 
                    Vol_II{II}=StackImage(cellfun(@(x) dicomread(x.Filename), InfosT1sIJ_II{II},'UniformOutput' , false)); 
                end
                aux=combnk([1:numel(Vol_II)],2);
                aux2=zeros(size(Vol_II));
                for II=1:size(aux,1)        % Remove repeated images
                    if isequal(Vol_II{aux(II,1)},Vol_II{aux(II,2)})
                        aux2(max(aux(II,1),aux(II,2)))=1;
                    end
                end
                InfosT1sIJ_II(logical(aux2))=[];
                %%%%%
                for II=1:numel(InfosT1sIJ_II)
                    Instances=cellfun(@(x) x.InstanceNumber, InfosT1sIJ_II{II}, 'UniformOutput', false);               
                    [~,si]=sort(cell2mat(Instances));
                    mkdir([WriteFolder PatientName filesep 'T1W' filesep Date filesep 'T1W_Alpha=' num2str(round(alpha))...
                        '_' num2str(J+II-1) filesep])
                    for K=1:length(si)
                        copyfile(InfosT1sIJ_II{II}{si(K)}.Filename,[WriteFolder PatientName filesep 'T1W' filesep Date...
                            filesep 'T1W_Alpha=' num2str(round(alpha)) '_' num2str(J+II-1) filesep num2str(K) '.dcm'])
                    end                    
                end
            end
        end
    end
end

% Visualize images, remove some if necesary

Images=AdjustDirVariable(dir([WriteFolder PatientName filesep 'T1W' filesep Date])); 
Images=Images([Images(:).isdir]);
Images=Images(cell2mat(cellfun(@(x) contains(x,'T1W_Alpha='), {Images(:).name}, 'UniformOutput', false)));
[Vols,Infs]=arrayfun(@(x) ReadDcmFolder4([x.folder filesep x.name filesep]), Images);


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
aux5=Images(kdisc);
if ~isempty(aux5)
    arrayfun(@(x) rmdir([x.folder filesep x.name filesep],'s'), aux5,'UniformOutput',false);    
end






function Out=ContrastPresent(Info)
Out=0;
try
    aux=Info.ContrastBolusAgent;
    if ~isempty(aux)
        Out=1;
    end
end



function KeyPressFcn (object, eventdata)

global Stop
keyPressed = eventdata.Key;
if strcmp(keyPressed,'escape')    
    Stop=1;
end




