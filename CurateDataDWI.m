function CurateDataDWI(WriteFolder,Infos,SeriesDescription,PatientName,Date,FolderDWI) 
% Separates DW data by different b-values        

clc
disp(['Curating DW Patient: ' PatientName ' Date: ' Date])
try
        rmdir([WriteFolder PatientName filesep 'DWI' filesep Date filesep],'s');
        rmdir([WriteFolder PatientName filesep 'DWI' filesep 'Registered' filesep Date filesep],'s');
end

aux2=cellfun(@(x) ~isempty(regexpi(x,FolderDWI)),...
        SeriesDescription, 'UniformOutput', false);      

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
    % Get b factor folder  
    try
        SNames=cellfun(@(x) x.SequenceName, InfosDWIsI, 'UniformOutput', false);
        [USNames,~,ISNames]=unique(SNames);    
        bs=cellfun(@(x) x(regexpi(x, '\d')), USNames, 'UniformOutput', false); 
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
    % For Each file identify b-factor
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
        aux=cellfun(@(x) strsplit(x,filesep), aux,'UniformOutput',false);
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
        mkdir([WriteFolder PatientName filesep 'DWI' filesep Date filesep 'b=' Bfact '_' num2str(I) filesep '0_' num2str(J)  filesep]);
        mkdir([WriteFolder PatientName filesep 'DWI' filesep Date filesep 'b=' Bfact '_' num2str(I) filesep 'Processed' filesep 'Registered' filesep '0_' num2str(J) filesep]);
        InfosJ=InfosDWIsI(b0s(J)==ids);
        for K=1:length(InfosJ)
            aux3=strsplit(InfosJ{K}.Filename,filesep);        
            copyfile(InfosJ{K}.Filename,[WriteFolder PatientName filesep 'DWI' filesep Date filesep 'b=' Bfact '_' num2str(I) filesep '0_' num2str(J)  filesep aux3{end}])
        end 
        copyfile([WriteFolder PatientName filesep 'DWI' filesep Date filesep 'b=' Bfact '_' num2str(I) filesep '0_' num2str(J)  filesep],...
            [WriteFolder PatientName filesep 'DWI' filesep Date filesep 'b=' Bfact '_' num2str(I) filesep 'Processed' filesep 'Registered' filesep '0_' num2str(J)  filesep])
    end
    for J=1:length(bN0s)
        mkdir([WriteFolder PatientName filesep 'DWI' filesep Date filesep 'b=' Bfact '_' num2str(I) filesep Bfact '_' num2str(J)  filesep]);
        mkdir([WriteFolder PatientName filesep 'DWI' filesep Date filesep 'b=' Bfact '_' num2str(I) filesep 'Processed' filesep 'Registered' filesep Bfact '_' num2str(J)  filesep]);
        InfosJ=InfosDWIsI(bN0s(J)==ids);
        for K=1:length(InfosJ)
            aux3=strsplit(InfosJ{K}.Filename,filesep);        
            copyfile(InfosJ{K}.Filename,[WriteFolder PatientName filesep 'DWI' filesep Date filesep 'b=' Bfact '_' num2str(I) filesep Bfact '_' num2str(J)  filesep aux3{end}])
        end   
        copyfile([WriteFolder PatientName filesep 'DWI' filesep Date filesep 'b=' Bfact '_' num2str(I) filesep Bfact '_' num2str(J)  filesep],...
            [WriteFolder PatientName filesep 'DWI' filesep Date filesep 'b=' Bfact '_' num2str(I) filesep 'Processed' filesep 'Registered' filesep Bfact '_' num2str(J)  filesep])
    end   
end

% Visualize images, remove some if necesary
Foldersb=AdjustDirVariable(dir([WriteFolder PatientName filesep 'DWI' filesep Date])); 
aux=arrayfun(@(x) AdjustDirVariable(dir([x.folder filesep x.name ])), Foldersb,'UniformOutput',false);    
aux2=vertcat(aux{:});    
aux3=arrayfun(@(x) strcmp(x.name,'Processed'), aux2);
aux2(aux3)=[];    
aux4=arrayfun(@(x) ReadDcmFolder4([x.folder filesep x.name filesep]), aux2);
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











    
   

