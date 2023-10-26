function MoveFilesForDateRegistrationT1(Folder,Date,ImageTemplate,Dates,ImageTemplates)

T1=[];AUC_NW=[];Ve_NW_NL=[];Vp_NW_NL=[];Kt_NW_NL=[];
if isequal(Date,Dates{1})    
    mkdir([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep])
    mkdir([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep])
    copyfile([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep ImageTemplate],...
        [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep ImageTemplate])
    copyfile([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep ImageTemplate],...
        [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep ImageTemplate])
    [~,Infos]=ReadDcmFolder4([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep ImageTemplate filesep]);
    load([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep 'T1Data_NoCorrection.mat'])
    try,load([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep '90sAUCs.mat']);end  % 9/15/2020
    try,load([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep 'Perfussion_Parameters_Maps.mat']);end
    Strs={'T1','AUC','Ve','Vp','Kt'}; 
    J=1;
    for PM={T1,AUC_NW,Ve_NW_NL,Vp_NW_NL,Kt_NW_NL}  
        if ~isempty(PM{1})
            PM{1}(imag(PM{1})>0)=0;
            WriteDicomFolderV2(real(PM{1}), Infos{1}, [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep Strs{J} filesep],Strs{J});   
            copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep Strs{J} filesep],...
                [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep Strs{J} filesep])
        end
        J=J+1;
    end


    MaskT1=load([Folder filesep 'T1W' filesep Date filesep 'Masks.mat']);  
    Labels=MaskT1.Labels;
    [VolT1,InfosT1]=ReadDcmFolder4([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep ImageTemplate filesep]);
    VolT1=VolT1{1};
    InfosT1=InfosT1{1};
    FieldsTumorT1=find(  (contains(Labels, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+[0-9]')), Labels)) ...
                        & ~contains(Labels, 'core', 'IgnoreCase',true) & ~contains(Labels, 'viable', 'IgnoreCase',true));
    MatMasksTumors=cat(4,MaskT1.MaksPerLabel{FieldsTumorT1});
    MatMasksTumors=any(MatMasksTumors,4);
    WriteDicomFolderV2(MatMasksTumors, InfosT1, [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'TumorSegmentationT1' filesep],...
        'TumorSegmentation',1,0) 
    copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'TumorSegmentationT1' filesep],...
        [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep 'TumorSegmentationT1' filesep])
else
     MaskT1_v1=load([Folder filesep 'T1W' filesep Dates{1} filesep 'Masks.mat']);  
     Labels_v1=MaskT1_v1.Labels;
     FieldsTumorT1_v1=find(  (contains(Labels_v1, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+[0-9]')), Labels_v1)) ...
                        & ~contains(Labels_v1, 'core', 'IgnoreCase',true) & ~contains(Labels_v1, 'viable', 'IgnoreCase',true));
                    
     MaskT1=load([Folder filesep 'T1W' filesep Date filesep 'Masks.mat']);  
     Labels=MaskT1.Labels;               
     FieldsTumorT1=find(  (contains(Labels, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+[0-9]')), Labels)) ...
                        & ~contains(Labels, 'core', 'IgnoreCase',true) & ~contains(Labels, 'viable', 'IgnoreCase',true));
                    
     CTT1_V1=regionprops(MaskT1_v1.MaksPerLabel{FieldsTumorT1_v1(1)},'Centroid');
     CTT1_V1=CTT1_V1.Centroid;
     CTT1=regionprops(MaskT1.MaksPerLabel{FieldsTumorT1(1)},'Centroid');
     CTT1=CTT1.Centroid;                     
                    
     [~,InfosT1_v1]=ReadDcmFolder4([Folder filesep 'T1W' filesep 'Registered' filesep Dates{1} filesep 'Local' filesep ImageTemplates{1} filesep]);
     InfosT1_v1=InfosT1_v1{1};  
     SizeT1_v1=size(MaskT1_v1.MaksPerLabel{FieldsTumorT1_v1(1)});
     SizeT1=size(MaskT1.MaksPerLabel{FieldsTumorT1(1)});
     [VolT1,InfosT1]=ReadDcmFolder4([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep ImageTemplate filesep]);
     InfosT1=InfosT1{1};
     VolT1=VolT1{1};
                    
     Delta=(CTT1_V1-((SizeT1_v1+1)/2)).*[InfosT1_v1{1}.PixelSpacing' InfosT1_v1{1}.SliceThickness];
     Xc_d=CTT1-Delta./[[InfosT1{1}.PixelSpacing(2) InfosT1{1}.PixelSpacing(1)]...
         InfosT1{1}.SliceThickness]; 
     Xd=((SizeT1+1)/2)-[Xc_d(2) Xc_d(1) Xc_d(3)];
     for I=1:3
        aux=[1:SizeT1(I)]+Xd(I);
        [aux2,Ind]=unique(round(aux));
        aux3=[1:SizeT1(I)];
        aux3=aux3(Ind);
        Ind=find(aux2>=1 & aux2<=SizeT1(I));
        aux2=aux2(Ind);
        aux3=aux3(Ind);
        IndsO{I}=aux3;  
        IndsD{I}=aux2;  
     end 
     Vol=zeros(SizeT1);
     Vol(IndsD{1},IndsD{2},IndsD{3})=VolT1(IndsO{1},IndsO{2},IndsO{3});
     
     mkdir([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep])
     mkdir([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep])
     WriteDicomFolderV2(Vol, InfosT1, [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep ImageTemplate filesep],...
        InfosT1{1}.SeriesDescription, 1,0)     
     
     copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep ImageTemplate],...
        [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep ImageTemplate])
    
    load([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep 'T1Data_NoCorrection.mat'])
    try,load([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep '90sAUCs.mat']);end  % 9/15/2020
    try,load([Folder filesep 'T1W' filesep 'Registered' filesep Date filesep 'Local' filesep 'Perfussion_Parameters_Maps.mat']);end
    Strs={'T1','AUC','Ve','Vp','Kt'}; 
    J=1;
    for PM={T1,AUC_NW,Ve_NW_NL,Vp_NW_NL,Kt_NW_NL}  
        if ~isempty(PM{1})
            Vol=zeros(SizeT1);
            Vol(IndsD{1},IndsD{2},IndsD{3})=PM{1}(IndsO{1},IndsO{2},IndsO{3});
            Vol(imag(Vol)>0)=0;
            WriteDicomFolderV2(real(Vol), InfosT1, [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep Strs{J} filesep],Strs{J});   
            copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep Strs{J} filesep],...
                [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep Strs{J} filesep])
        end
        J=J+1;
    end
    MatMasksTumors=cat(4,MaskT1.MaksPerLabel{FieldsTumorT1});
    MatMasksTumors=any(MatMasksTumors,4);
    Vol=zeros(SizeT1);
    Vol(IndsD{1},IndsD{2},IndsD{3})=MatMasksTumors(IndsO{1},IndsO{2},IndsO{3});
    WriteDicomFolderV2(Vol, InfosT1, [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'TumorSegmentationT1' filesep],...
        'TumorSegmentation',1,0) 
    copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'TumorSegmentationT1' filesep],...
        [Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep 'TumorSegmentationT1' filesep])  
     
    
end








