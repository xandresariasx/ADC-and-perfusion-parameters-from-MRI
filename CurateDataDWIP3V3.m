
function CurateDataDWIP3V3(WriteFolder,Infos,SeriesDescription,PatientName,Date,FolderDWI,BiasCorrect) 

clc
disp(['Curating DW Patient: ' PatientName ' Date: ' Date])
try
rmdir([WriteFolder PatientName '\DWI\' Date '\'],'s');
rmdir([WriteFolder PatientName '\DWI\Registered\' Date '\'],'s');
end

aux2=cellfun(@(x) ~isempty(regexpi(x,FolderDWI)),...
        SeriesDescription, 'UniformOutput', false);      % For SMIL only

aux2=cell2mat(aux2);

InfosDWI=Infos(aux2);
SeriesDescriptionDWI=SeriesDescription(aux2);

if isempty(InfosDWI)
    return
end

% Get number of folders in Date
IDs=cellfun(@(x) x.SeriesInstanceUID, InfosDWI, 'UniformOutput', false);  
[Ut1,~,It1]=unique(IDs);
 
 
for I=1:length(Ut1)
    InfosDWIsI=InfosDWI(I==It1);
    SeriesDescriptionDWIsI=SeriesDescriptionDWI(I==It1);
    
    % Find Number of Images in Folder
    Instances=cellfun(@(x) x.InstanceNumber, InfosDWIsI, 'UniformOutput', false);
    SL=cell2mat(cellfun(@(x) x.SliceLocation, InfosDWIsI, 'UniformOutput', false)); 
    NImages=numel(Instances)/numel(unique(SL));
    
    % Get b factor folder  (Check!!!)
    try
        SNames=cellfun(@(x) x.SequenceName, InfosDWIsI, 'UniformOutput', false);
        [USNames,~,ISNames]=unique(SNames);    
        bs=cellfun(@(x) x(regexpi(x,'\d')), USNames, 'UniformOutput', false); 
        if isempty(bs{1})
            GenerateError=1/bs;   % Force error in case bs is empty
        end
        bs=sort(bs);
        Bfact=bs{end};  
    catch
        bs=cellfun(@(x) x(regexpi(x,'\d')), SeriesDescriptionDWIsI, 'UniformOutput', false); 
        bs=cellfun(@(x) x(1:min([numel(x) 3])), bs, 'UniformOutput', false);
        Bfact=unique(bs);  
        Bfact=Bfact{end};   
        if numel(Bfact)<3
            Bfact='450';
        end
    end   
    
    % for Each file identify b-factor
    if numel(unique(cell2mat(Instances)))~=numel(cell2mat(Instances))
        aux1=unique(cell2mat(Instances));
        Vol1=StackImage(cellfun(@(x) dicomread(x.Filename), InfosDWIsI(1:max(aux1)),'UniformOutput' , false));  
        Vol2=StackImage(cellfun(@(x) dicomread(x.Filename), InfosDWIsI(max(aux1)+1:end),'UniformOutput' , false));
        if isequal(Vol1,Vol2)
            InfosDWIsI(max(aux1)+1:end)=[];
            Instances(max(aux1)+1:end)=[];
        end
    end
    if numel(unique(cell2mat(Instances)))~=numel(cell2mat(Instances))
        aux=cellfun(@(x) x.Filename, InfosDWIsI,'UniformOutput',false);
        aux=cellfun(@(x) strsplit(x,'\'), aux,'UniformOutput',false);
        aux=cellfun(@(x) strsplit(x{end},'-'), aux,'UniformOutput',false);
        aux=cellfun(@(x) x{2}, aux, 'UniformOutput',false);
        Vols=[];
        ids=[];
        for ID=unique(aux)
            ID=ID{1};
            InfosID=InfosDWIsI(strcmp(aux,ID));
            InstancesID=cellfun(@(x) x.InstanceNumber, InfosID, 'UniformOutput', false);
            [~,si]=sort(cell2mat(InstancesID));
            InfosIDsI=InfosID(si);
            Vol=StackImage(cellfun(@dicomread, InfosIDsI, 'UniformOutput', false));
            N=numel(InfosIDsI);
            SL=cell2mat(cellfun(@(x) x.SliceLocation, InfosIDsI, 'UniformOutput', false)); 
            NImages=numel(InstancesID)/numel(unique(SL));
            if isempty(ids)
                Mid=0;Lid=0;
            else
                Mid=max(ids);
                Lid=numel(ids);
            end        
            for J=1:NImages, 
                Vols{end+1}=Vol(:,:,1+(N/NImages)*(J-1):J*N/NImages); 
                ids(Lid+1+(N/NImages)*(J-1):Lid+J*N/NImages)= J+Mid;   
            end          
        end 
        Means=cellfun(@(x) mean(x(:)), Vols);
        b0s=find(arrayfun(@(x) any(x>1.5*Means), Means));
        bN0s=find(~arrayfun(@(x) any(x>1.5*Means), Means)); 
    else
        [~,si]=sort(cell2mat(Instances));
        InfosDWIsI=InfosDWIsI(si);
        Vol=StackImage(cellfun(@dicomread, InfosDWIsI, 'UniformOutput', false));
        N=numel(InfosDWIsI);
        ids=zeros(size(InfosDWIsI));
        SL=cell2mat(cellfun(@(x) x.SliceLocation, InfosDWIsI, 'UniformOutput', false)); 
        NImages=numel(Instances)/numel(unique(SL));
        Vols=[];
        for J=1:NImages, 
            Vols{J}=Vol(:,:,1+(N/NImages)*(J-1):J*N/NImages); 
            ids(1+(N/NImages)*(J-1):J*N/NImages)= J;   
        end
        Means=cellfun(@(x) mean(x(:)), Vols);
        b0s=find(arrayfun(@(x) any(x>1.5*Means), Means));
        bN0s=find(~arrayfun(@(x) any(x>1.5*Means), Means));
    end
    
    % Create folders curated images
    for J=1:length(b0s)
        mkdir([WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\0_' num2str(J)  '\']);
        mkdir([WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\Processed\Registered\0_' num2str(J) '\']);
        InfosJ=InfosDWIsI(b0s(J)==ids);
        for K=1:length(InfosJ)
            aux3=strsplit(InfosJ{K}.Filename,'\');        
            copyfile(InfosJ{K}.Filename,[WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\0_' num2str(J)  '\' aux3{end}])
        end    
        if BiasCorrect
            system(['C:\Temp\N4\bin\Release\N4 '...
                [WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\0_' num2str(J)] ' temp.mha'])
            mhaToDicom([WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\0_' num2str(J)  '\']);                      
        end
        copyfile([WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\0_' num2str(J)  '\'],...
            [WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\Processed\Registered\0_' num2str(J)  '\'])
    end
    for J=1:length(bN0s)
        mkdir([WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\' Bfact '_' num2str(J)  '\']);
        mkdir([WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\Processed\Registered\' Bfact '_' num2str(J)  '\']);
        InfosJ=InfosDWIsI(bN0s(J)==ids);
        for K=1:length(InfosJ)
            aux3=strsplit(InfosJ{K}.Filename,'\');        
            copyfile(InfosJ{K}.Filename,[WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\' Bfact '_' num2str(J)  '\' aux3{end}])
        end    
        if BiasCorrect
            system(['C:\Temp\N4\bin\Release\N4 '...
                [WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\' Bfact '_' num2str(J)] ' temp.mha'])
            mhaToDicom([WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\' Bfact '_' num2str(J)  '\']);
        end
        copyfile([WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\' Bfact '_' num2str(J)  '\'],...
            [WriteFolder PatientName '\DWI\' Date '\b=' Bfact '_' num2str(I) '\Processed\Registered\' Bfact '_' num2str(J)  '\'])
    end   
end



% Visualize images, remove some if necesary

Foldersb=AdjustDirVariable(dir([WriteFolder PatientName '\DWI\' Date])); 
aux=arrayfun(@(x) AdjustDirVariable(dir([x.folder '\' x.name ])), Foldersb,'UniformOutput',false);    
aux2=vertcat(aux{:});    
aux3=arrayfun(@(x) strcmp(x.name,'Processed'), aux2);
aux2(aux3)=[];    
aux4=arrayfun(@(x) ReadDcmFolder4([x.folder '\' x.name '\']), aux2);
aux6=cell2mat(aux4);
close all force
Fig=figure('units','normalized','outerposition',[0 0 1 1]); colormap('gray')
for I=1:length(aux4)
    subplot(ceil(length(aux4)/4),4,I), imagesc(aux4{I}(:,:,floor(size(aux4{I},3)/2)),...
        [0 prctile(aux6(:),99)]) 
    title([aux2(I).name '   (' num2str(I) ')'],'Interpreter','none')
    axis equal
    axis off
end
set (Fig, 'WindowKeyPressFcn', @KeyPressFcn);
cont=1;
global Stop
Stop=0;
while Stop==0 
    for I=1:length(aux4)
        if cont<=size(aux4{I},3)
            subplot(ceil(length(aux4)/4),4,I), imagesc(aux4{I}(:,:,cont),...
                [0 prctile(aux6(:),99)])     
            title([aux2(I).name '   (' num2str(I) ')' '  (Slice #: ' num2str(cont) ' )'],'Interpreter','none')
            axis equal
            axis off    
        end
    end
    pause(0.5)
    cont=cont+1;
    if cont>max(cellfun(@(x) size(x,3), aux4))
        cont=1;
    end
end
for I=1:length(aux4)
    subplot(ceil(length(aux4)/4),4,I), imagesc(aux4{I}(:,:,floor(size(aux4{I},3)/2)),...
        [0 prctile(aux6(:),99)]) 
    title([aux2(I).name '   (' num2str(I) ')'],'Interpreter','none')
    axis equal
    axis off
end

aux=combnk([1:numel(aux4)],2);
RM=[];
for I=1:size(aux,1)
    if isequal(aux4{aux(I,1)},aux4{aux(I,2)})
        RM=[RM max(aux(I,1),aux(I,2))];
    end
end
if ~isempty(RM)
    disp(['Automatically discarded volumes: ' num2str(RM)])
end
kdisc=str2num(input('Discard Other Volumes? ','s'));
kdisc=unique([kdisc RM]);
aux5=aux2(kdisc);


for I=1:length(aux4)
    subplot(ceil(length(aux4)/4),4,I), imagesc(aux4{I}(:,:,floor(size(aux4{I},3)/2)),...
        [0 prctile(aux6(:),99)]) 
    if any(kdisc==I)
        title([aux2(I).name '   (' num2str(I) ', REMOVED)'],'Interpreter','none')
    else
        title([aux2(I).name '   (' num2str(I) ')'],'Interpreter','none')
    end
    axis equal
    axis off
end
savefig([WriteFolder PatientName '\DWI\' Date '\DWImages'])
saveas(gcf,[WriteFolder PatientName '\DWI\' Date '\DWImages.jpeg'])


if ~isempty(aux5)
    arrayfun(@(x) rmdir([x.folder '\' x.name '\'],'s'), aux5,'UniformOutput',false);
    arrayfun(@(x) rmdir([x.folder '\Processed\Registered\' x.name '\'],'s'), aux5,'UniformOutput',false);  
    Foldersb=AdjustDirVariable(dir([WriteFolder PatientName '\DWI\' Date])); 
    aux=arrayfun(@(x) AdjustDirVariable(dir([x.folder '\' x.name ])), Foldersb,'UniformOutput',false);    
    aux2=vertcat(aux{:}); 
    [aux7,~,aux9]=unique({aux2(:).folder});
    aux10=hist(aux9,[1:max(aux9)]);
    aux11=aux7(aux10==1);
    if ~isempty(aux11)
        cellfun(@(x) rmdir(x,'s'), aux11,'UniformOutput',false);
    end
end

fid=fopen([WriteFolder PatientName '\DWI\' Date '\DiscardedDW.txt'],'w');
fprintf(fid, num2str(kdisc));
fclose(fid);




function KeyPressFcn (object, eventdata)

global Stop
keyPressed = eventdata.Key;
if strcmp(keyPressed,'escape')
    aux=get(gca,'children');
    delete(aux(1))
    Stop=1;
end









    
   

