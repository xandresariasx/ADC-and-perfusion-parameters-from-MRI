function MoveFilesForDWT1RegistrationV2(Folder,Dates,LocalRegTemplateB,ImageTemplateB,RotateDW)
            
            
try,rmdir([Folder filesep 'T1-DW Registered' filesep],'s'),end
mkdir([Folder filesep 'T1-DW Registered' filesep])
mkdir([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep])
mkdir([Folder filesep 'T1-DW Registered' filesep 'T1W' filesep])
aux=AdjustDirVariable(dir([Folder filesep 'DWI' filesep...
    'DateRegisteredLocal' filesep Dates{1} filesep 'Processed' filesep]));
% Move template DW from baseline
copyfile([Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Dates{1} filesep 'Processed' filesep aux(end).name filesep LocalRegTemplateB],...
    [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep LocalRegTemplateB])
% Move ADC from local date registration
for Date=Dates
    copyfile([Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Date{1} filesep 'Processed' filesep aux(end).name filesep 'ADC' filesep],...
        [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'ADC_' Date{1}])
end
% save tumor segmentations in DW baseline as an image
copyfile([Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Dates{1} filesep 'Processed' filesep aux(end).name filesep 'TumorSegmentationDW' filesep],...
    [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'TumorSegmentationDW' filesep])

% MaskDW=load([Folder filesep 'DWI' filesep 'Registered' filesep Dates{1} filesep 'Global' filesep 'Masks.mat']);  
% Labels=MaskDW.Labels;
[VolDW,InfosDW]=ReadDcmFolder4([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep LocalRegTemplateB filesep]);
VolDW=VolDW{1};
InfosDW=InfosDW{1};
% FieldsTumorDW=find(  (contains(Labels, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+[0-9]')), Labels)) ...
%                     & ~contains(Labels, 'core', 'IgnoreCase',true) & ~contains(Labels, 'viable', 'IgnoreCase',true));
% MatMasksTumors=cat(4,MaskDW.MaksPerLabel{FieldsTumorDW});
% MatMasksTumors=any(MatMasksTumors,4);
% WriteDicomFolderV2(MatMasksTumors, InfosDW, [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'TumorSegmentationDW' filesep],...
%     'TumorSegmentation',1,0) 

if ~isempty(RotateDW)
    if strcmp(RotateDW,'AxToSag')
        z=RotateAxToSag([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep]);
        [VolDW,InfosDW]=ReadDcmFolder4([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep LocalRegTemplateB filesep]);
        VolDW=VolDW{1};
        InfosDW=InfosDW{1};
    end
     if strcmp(RotateDW,'AxToCor')
        z=RotateAxToCor([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep]);
        [VolDW,InfosDW]=ReadDcmFolder4([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep LocalRegTemplateB filesep]);
        VolDW=VolDW{1};
        InfosDW=InfosDW{1};
    end
end

% Move template baseline T1
copyfile([Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Dates{1} filesep ImageTemplateB],...
    [Folder filesep 'T1-DW Registered' filesep 'T1W' filesep ImageTemplateB])
copyfile([Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Dates{1} filesep ImageTemplateB],...
    [Folder filesep 'T1-DW Registered' filesep 'T1W' filesep 'Processed' filesep 'Registered' filesep ImageTemplateB])

% Shift DW images to be at similar location to T1
MaskT1=load([Folder filesep 'T1W' filesep Dates{1} filesep 'Masks.mat']);  % Use the location of the first tumor as reference
MaskDW=load([Folder filesep 'DWI' filesep 'Registered' filesep Dates{1} filesep 'Global' filesep 'Masks.mat']); 
FieldsTumorDW=find(  (contains(MaskDW.Labels, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+1')), MaskDW.Labels)) ...
                    & ~contains(MaskDW.Labels, 'core', 'IgnoreCase',true) & ~contains(MaskDW.Labels, 'viable', 'IgnoreCase',true)); 
FieldsTumorT1=find(  (contains(MaskT1.Labels, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+1')), MaskT1.Labels)) ...
                    & ~contains(MaskT1.Labels, 'core', 'IgnoreCase',true) & ~contains(MaskT1.Labels, 'viable', 'IgnoreCase',true));    

if ~isempty(RotateDW)
    if strcmp(RotateDW,'AxToSag')
        MaskDW.MaksPerLabel{FieldsTumorDW(1)}=permute(MaskDW.MaksPerLabel{FieldsTumorDW(1)},[3 1 2]);
        MaskDW.MaksPerLabel{FieldsTumorDW(1)}=MaskDW.MaksPerLabel{FieldsTumorDW(1)}...
            (end:-1:1,:,z:size(MaskDW.MaksPerLabel{1},3)+1-z);
    end
    if strcmp(RotateDW,'AxToCor')
        MaskDW.MaksPerLabel{FieldsTumorDW(1)}=permute(MaskDW.MaksPerLabel{FieldsTumorDW(1)},[3 2 1]);
        MaskDW.MaksPerLabel{FieldsTumorDW(1)}=MaskDW.MaksPerLabel{FieldsTumorDW(1)}...
            (end:-1:1,:,z:size(MaskDW.MaksPerLabel{1},3)+1-z);
    end
end                
CTDW=regionprops(MaskDW.MaksPerLabel{FieldsTumorDW(1)},'Centroid');
CTDW=CTDW.Centroid;
CTT1=regionprops(MaskT1.MaksPerLabel{FieldsTumorT1(1)},'Centroid');
CTT1=CTT1.Centroid;               
[~,InfosT1]=ReadDcmFolder4([Folder filesep 'T1-DW Registered' filesep 'T1W' filesep ImageTemplateB filesep]);
InfosT1=InfosT1{1};
SizeT1=size(MaskT1.MaksPerLabel{FieldsTumorT1(1)});
sizeDW=size(VolDW);
% Delta=(CTT1(3)-((SizeT1(3)+1)/2))*InfosT1{1}.SliceThickness;  % Number of slices from tumor center to image center in T1
% Xc_d=CTDW(3)-Delta/InfosDW{1}.SliceThickness;  % Slice location with same anatomical position as T1 image center
% zd=((sizeDW(3)+1)/2)-Xc_d;  % Vertical shift
% aux=[1:sizeDW(3)]+zd;  % Vertical shift transformation
% [aux2,Ind]=unique(round(aux));
% aux3=[1:sizeDW(3)];
% aux3=aux3(Ind);
% Ind=find(aux2>=1 & aux2<=sizeDW(3));
% aux2=aux2(Ind);
% aux3=aux3(Ind);
% Vol=zeros(size(VolDW));
% Vol(:,:,aux2)=VolDW(:,:,aux3);

Delta=(CTT1-((SizeT1+1)/2)).*[InfosT1{1}.PixelSpacing' InfosT1{1}.SliceThickness];
Xc_d=CTDW-Delta./[[InfosDW{1}.PixelSpacing(2) InfosDW{1}.PixelSpacing(1)]...
    InfosDW{1}.SliceThickness]; 
Xd=((sizeDW+1)/2)-[Xc_d(2) Xc_d(1) Xc_d(3)];
for I=1:3
    aux=[1:sizeDW(I)]+Xd(I);
    [aux2,Ind]=unique(round(aux));
    aux3=[1:sizeDW(I)];
    aux3=aux3(Ind);
    Ind=find(aux2>=1 & aux2<=sizeDW(I));
    aux2=aux2(Ind);
    aux3=aux3(Ind);
    IndsO{I}=aux3;  
    IndsD{I}=aux2;  
end
Vol=zeros(size(VolDW));
Vol(IndsD{1},IndsD{2},IndsD{3})=VolDW(IndsO{1},IndsO{2},IndsO{3});

rmdir([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep LocalRegTemplateB filesep],'s')
WriteDicomFolderV2(Vol, InfosDW, [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep LocalRegTemplateB filesep],...
    InfosDW{1}.SeriesDescription, 1,0)
copyfile([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep LocalRegTemplateB filesep],...
    [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'Processed' filesep 'Registered' filesep LocalRegTemplateB filesep])

[Vol,Infos]=ReadDcmFolder4([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'TumorSegmentationDW' filesep]);
Vol=Vol{1};Infos=Infos{1}; 
rmdir([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'TumorSegmentationDW' filesep],'s')
Vol2=zeros(size(VolDW));
Vol2(IndsD{1},IndsD{2},IndsD{3})=Vol(IndsO{1},IndsO{2},IndsO{3});
WriteDicomFolderV2(Vol2, Infos, [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'TumorSegmentationDW' filesep],...
    Infos{1}.SeriesDescription, 1,0)
copyfile([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'TumorSegmentationDW' filesep],...
    [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'Processed' filesep 'Registered' filesep 'TumorSegmentationDW' filesep])
for Date=Dates
    [Vol,Infos]=ReadDcmFolder4([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'ADC_' Date{1} filesep]);
    Vol=Vol{1};Infos=Infos{1};   
    rmdir([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'ADC_' Date{1} filesep],'s')
    Vol2=zeros(size(VolDW));
    Vol2(IndsD{1},IndsD{2},IndsD{3})=Vol(IndsO{1},IndsO{2},IndsO{3});
    WriteDicomFolderV2(Vol2, Infos, [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'ADC_' Date{1} filesep],...
        Infos{1}.SeriesDescription)
    copyfile([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'ADC_' Date{1} filesep],...
        [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'Processed' filesep 'Registered' filesep 'ADC_' Date{1} filesep])
end


function zc=RotateAxToSag(Folder)

Imgs=AdjustDirVariable(dir(Folder) ); 

[VolOrg,FilesOrg]=ReadDcmFolder3([Folder Imgs(end-1).name filesep]);
VolOrg=VolOrg{1};
VolOrgSag=permute(VolOrg,[3 1 2]);
VolOrgSag=VolOrgSag(end:-1:1,:,:);
% [~,~,zc]=ind2sub(size(VolOrgSag),find(VolOrgSag,1));
zc=1;

for I=1:numel(Imgs)    
    SrcFolder=[Folder Imgs(I).name filesep];    
    [VolOrg,FilesOrg]=ReadDcmFolder3(SrcFolder);
    VolOrg=VolOrg{1};
    VolOrgSag=permute(VolOrg,[3 1 2]);
    VolOrgSag=VolOrgSag(end:-1:1,:,zc:size(VolOrgSag,3)-zc);

    zr=FilesOrg{1}{1}.SpacingBetweenSlices;
    xryr=FilesOrg{1}{1}.PixelSpacing;
    xr=xryr(1);
    yr=xryr(2);
    z=yr;
    x=zr;
    y=xr;
    rmdir(SrcFolder,'s')
    DestImage=SrcFolder;
    mkdir(DestImage)
    for I=1:size(VolOrgSag,3)
        dicomwrite(uint32(VolOrgSag(:,:,I)), [DestImage ...
            num2str(I) '.dcm'],   'MultiframeSingleFile',false)    
        info=dicominfo([DestImage num2str(I) '.dcm']);
        info.SpacingBetweenSlices=z;
        info.SliceLocation=(I-1)*z;
        info.SliceThickness=z;
        info.PixelSpacing=[x;y];
        info.InstanceNumber=I;
        info.SeriesInstanceUID=FilesOrg{1}{1}.SeriesInstanceUID;
        info.SeriesDescription=FilesOrg{1}{1}.SeriesDescription; 
        try
            info.RescaleIntercept=FilesOrg{1}{1}.RescaleIntercept;           
            info.RescaleSlope=FilesOrg{1}{1}.RescaleSlope;   
        end
        dicomwrite(uint32(VolOrgSag(:,:,I)), [DestImage ...
            num2str(I) '.dcm'], info, 'MultiframeSingleFile',false,'CreateMode','copy')
    end
end


function zc=RotateAxToCor(Folder)

Imgs=AdjustDirVariable(dir(Folder) ); 

[VolOrg,FilesOrg]=ReadDcmFolder3([Folder Imgs(end-1).name filesep]);
VolOrg=VolOrg{1};
VolOrgCor=permute(VolOrg,[3 2 1]);
VolOrgCor=VolOrgCor(end:-1:1,:,:);
% [~,~,zc]=ind2sub(size(VolOrgSag),find(VolOrgSag,1));
% zc=zc-2;
zc=1;

for I=1:numel(Imgs)    
    SrcFolder=[Folder Imgs(I).name filesep];    
    [VolOrg,FilesOrg]=ReadDcmFolder3(SrcFolder);
    VolOrg=VolOrg{1};
    VolOrgCor=permute(VolOrg,[3 2 1]);
    VolOrgCor=VolOrgCor(end:-1:1,:,:);
    VolOrgCor=VolOrgCor(:,:,zc:size(VolOrgCor,3)+1-zc);

    zr=FilesOrg{1}{1}.SpacingBetweenSlices;
    xryr=FilesOrg{1}{1}.PixelSpacing;
    xr=xryr(1);
    yr=xryr(2);
    z=xr;
    x=zr;
    y=yr;
    rmdir(SrcFolder,'s')
    DestImage=SrcFolder;
    mkdir(DestImage)
    for I=1:size(VolOrgCor,3)
        dicomwrite(uint32(VolOrgCor(:,:,I)), [DestImage ...
            num2str(I) '.dcm'],   'MultiframeSingleFile',false)    
        info=dicominfo([DestImage num2str(I) '.dcm']);
        info.SpacingBetweenSlices=z;
        info.SliceLocation=(I-1)*z;
        info.SliceThickness=z;
        info.PixelSpacing=[x;y];
        info.InstanceNumber=I;
        info.SeriesInstanceUID=FilesOrg{1}{1}.SeriesInstanceUID;
        info.SeriesDescription=FilesOrg{1}{1}.SeriesDescription; 
        try
            info.RescaleIntercept=FilesOrg{1}{1}.RescaleIntercept;           
            info.RescaleSlope=FilesOrg{1}{1}.RescaleSlope;   
        end
        dicomwrite(uint32(VolOrgCor(:,:,I)), [DestImage ...
            num2str(I) '.dcm'], info, 'MultiframeSingleFile',false,'CreateMode','copy')
    end
end



















            
            
            