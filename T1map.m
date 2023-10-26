function T1map(Folder,FolderMask,RefTiss, ScaleT1Ims, gammaT1, gammaM0,...
    CorrectBiasBetweenSlices, AssumeT1s, AirTh, solver,Date,Dates,WriteFolder,...
    PatientName,CorrectFA,WeightT1)

if CorrectFA
    T1mapCorrectFA(Folder,FolderMask,gammaT1, gammaM0,...
        solver,AirTh,Date,Dates,WriteFolder,PatientName,WeightT1)
    return
end

VolAvg=GetAvDCEImg(FolderMask);
Air=VolAvg/max(VolAvg(:))<AirTh;   % Air in the image

Dates=Dates(cellfun(@(x) exist([WriteFolder PatientName filesep 'T1W' filesep x filesep])~=0,Dates));

if ScaleT1Ims
    if 1%strcmp(Date,Dates{1})
        [Vols,Info, C30]=ReScaleT1ImagesV4(Folder, FolderMask, RefTiss);
    else
        try                                         % To use the same scaling factors as visit 1
            [Vols,Info]=ReadVolsT1(Folder);
        catch
            CorrectTagsDicom(Folder);
            [Vols,Info]=ReadVolsT1(Folder);
        end
        SD=load([WriteFolder PatientName filesep 'T1W' filesep 'Registered' filesep Dates{1} filesep 'Local' filesep 'ScalingData2.mat']); %,'Cs','AnglesD','trS')         
        trS=cellfun(@(x) x{1}.RepetitionTime, Info)/1000;  
        AnglesD=cellfun(@(x) x{1}.FlipAngle, Info);
        [~,Order]=sort(AnglesD);
        AnglesD=AnglesD(Order);
        trS=trS(Order);   
        Info=Info(Order);
        Vols=Vols(Order);
        if isequal(AnglesD,SD.AnglesD) && isequal(trS,SD.trS)
            Vols=cellfun(@(x,y) x./y, Vols, num2cell(SD.Cs),'UniformOutput',false);
            C30=SD.Cs;
        else
            [Vols,Info, C30]=ReScaleT1ImagesV4(Folder, FolderMask, RefTiss);
        end   
    end
else
    try
        [Vols,Info]=ReadVolsT1(Folder);
    catch
        CorrectTagsDicom(Folder);
        [Vols,Info]=ReadVolsT1(Folder);
    end
    if CorrectBiasBetweenSlices
        Vols=cellfun(@(x) CorrectBiasBetweenSlicesF(x), Vols,'UniformOutput',false);
    end
    C30=1;
    ErrorOri=[]; val=[]; Method='Original';  RefTiss=[];
    save([Folder 'ScalingData.mat'],'ErrorOri','val','Method','RefTiss')
end
Air=Air | ~all(cat(4,Vols{:}),4);  % 11/30/2020 if T1w is undefined is also air

%%%
% Overlay3D(VolAvg,...
%     Air,0.3,[0 1],0,[])
%%%

warning off

AnglesR=[];
try
trS=cellfun(@(x) x{1}.RepetitionTime, Info)/1000;   % from ms to s
tr=mean(trS);
AnglesD=cellfun(@(x) x{1}.FlipAngle, Info);
AnglesR=AnglesD*pi/180;                                         % from  degrees to radians
end

%%% Assume T1 and Mo if no T1W images or only 1 FA
if isempty(AnglesR) | numel(unique(AnglesR))<=1 | AssumeT1s
    if ~isempty(Info) && ~AssumeT1s
        T1=AssumeT1([FolderMask 'Masks.mat'],Info{1}{1}.MagneticFieldStrength);
        A=((1-exp(-trS(1)/T1))./(1-cos(AnglesR(1)).*exp(-trS(1)/T1))).*sin(AnglesR(1));
        %M0=mean(cat(4,Vols{:}),4)./ A;
        M0=Vols{1}./ A;
    else
        T1=AssumeT1([FolderMask 'Masks.mat'],1.5);
        [Vols,Info]=cellfun(@(x) ReadDcmFolder4([Folder x filesep]),...
            {'DCE_t=1','DCE_t=2','DCE_t=3'});
        Vols=cat(4,Vols{:});
        Vol=mean(Vols,4); 
        trS=cellfun(@(x) x{1}.RepetitionTime, Info)/1000;   % from ms to s
        AnglesD=cellfun(@(x) x{1}.FlipAngle, Info);
        AnglesR=AnglesD*pi/180; 
        A=((1-exp(-trS(1)/T1))./(1-cos(AnglesR(1)).*exp(-trS(1)/T1))).*sin(AnglesR(1));
        M0=Vol./ A;
    end
    
    T1(Air)=0;
    M0(Air)=0;
    save([Folder filesep 'T1Data_NoCorrection.mat'],'T1','M0')
    save([Folder filesep 'C30.mat'],'C30')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gammaT1 = 5e-5;
% gammaM0 = 5e-5;
h=[Info{1}{1}.PixelSpacing' Info{1}{1}.SpacingBetweenSlices];
%%% Crop Volumes for speed
Stats=regionprops3(~Air); 
aux=Stats.BoundingBox;
aux(:,4:end)=aux(:,4:end)+aux(:,1:3)-1;
if size(aux,1)==1
    xi=round((aux(:,1:3)));
    xl=round((aux(:,4:end)));
else
    xi=round(min(aux(:,1:3)));
    xl=round(max(aux(:,4:end)));
end
[xi,xl]=CorrectXiXl(xi,xl,Air);   % 1/28/2021 To guarantee size(voi)>3
VolsC=cellfun(@(x) x(xi(2):xl(2),xi(1):xl(1),xi(3):xl(3)),Vols,'UniformOutput',false);
AirC=Air(xi(2):xl(2),xi(1):xl(1),xi(3):xl(3));
%%%
[T1c,M0c] = getNLParameterfitDiffV2(cat(4,VolsC{:}),h,AnglesR,trS*1000,...
    gammaT1,gammaM0,AirC,[],[],WeightT1,'solver',solver);
T1=zeros(size(Air));
T1(xi(2):xl(2),xi(1):xl(1),xi(3):xl(3))=T1c;
M0=zeros(size(Air));
M0(xi(2):xl(2),xi(1):xl(1),xi(3):xl(3))=M0c;
T1(Air)=0;        % 9/22/2020
M0(Air)=0;
T1=T1/1000;


save([Folder filesep 'T1Data_NoCorrection.mat'],'T1','M0')

save([Folder filesep 'C30.mat'],'C30')

figure('units','normalized','outerposition',[0 0 1 1]), colormap('gray')
Maps={T1, M0};
Limits={[0 3], [0 10000], [0 3]};
Names={'T1', 'Mo'};
load([FolderMask 'Masks.mat'])
SlicePos=GetSlicePos(Labels,MaksPerLabel,size(T1));
for I=1:length(Maps)
    subplot(1,2,I), 
    Overlay(Vols{1}(:,:,SlicePos), Maps{I}(:,:,SlicePos), 0.5, Limits{I})    
    title(Names{I},'Interpreter','none')
end
savefig([Folder filesep 'T1ParMaps_NoCorrection'])
saveas(gcf,[Folder filesep 'T1ParMaps_NoCorrection.jpeg'])



function Vol=GetAvDCEImg(Folder)

Files=dir([Folder filesep 'Processed' filesep 'Registered' filesep]);
FileNames={Files(:).name};
aux=cellfun(@(x) contains(x,'DCE_t='), FileNames);
if all(aux==0)
    aux=cellfun(@(x) contains(x,'Alpha='), FileNames);
end
FileNames=FileNames(aux);
Vols=cellfun(@(x) ReadDcmFolder4([Folder filesep 'Processed' filesep 'Registered' filesep x filesep]),...
    FileNames);
Vols=cat(4,Vols{:});
Vol=mean(Vols,4);  


function [VolsRT1,Info]=ReadVolsT1(FolderT1Vols)

Files=AdjustDirVariable(dir([FolderT1Vols 'T1W_Alpha=*']));
if isempty(Files)
    VolsRT1=[];Info=[];
    return
end
VolsRT1=[];Info=[];
for I=1:length(Files)
    [VolsRT1{I},Info{I}]=ReadDcmFolder3([FolderT1Vols Files(I).name filesep]);
    VolsRT1{I}=VolsRT1{I}{1}; Info{I}=Info{I}{1};
end
trS=cellfun(@(x) x{1}.RepetitionTime, Info)/1000;  
AnglesD=cellfun(@(x) x{1}.FlipAngle, Info);
[~,Order]=sort(AnglesD);
AnglesD=AnglesD(Order);
trS=trS(Order);
Info=Info(Order);
VolsRT1=VolsRT1(Order);


function T1=AssumeT1(FileMask,FS)

load(FileMask)
TableT1L={'Tissue',                      '1.5T',       '3T';...
           'tumor',                         0.7,        0.8;...
           'Muscle',                          1,        1.4;...
           'liver',                        0.57,       0.8;...
           'kidney/RenalCortex',            0.9,        1.1;...
           'fat',                          0.34,       0.35;...
           'bone/BoneMarrow',              0.54,        0.58;...
           'spleen',                          1,        1.3;...
           'aorta/vein/blood/artery',       1.4,        1.9};
       
T1=ones(size(MaksPerLabel{1}));
for I=1:numel(Labels)
    TI=find(contains(TableT1L(2:end,1), Labels{I},'IgnoreCase',true),1)+1;
    FI=find(strcmp(TableT1L(1,2:3),[num2str(FS) 'T']),1)+1;
    T1v=TableT1L{TI,FI};
    T1(MaksPerLabel{I}==1)=T1v;
end


function Vol_c=CorrectBiasBetweenSlicesF(Vol)

aux=Vol(:,:,floor(size(Vol,3)/2));
aux=aux(:);
Is95=prctile(aux,95);

for I=1:size(Vol,3)
    aux=Vol(:,:,I);
    I95=prctile(aux(:),95);
    Vol_c(:,:,I)=(Is95/I95)*Vol(:,:,I);
end

function [xi,xl]=CorrectXiXl(xi,xl,Air)

 for I=1:3
    if abs(xi(I)-xl(I))<3
        if xi(I)-1>0
            xi(I)=xi(I)-1;
        else
            if xl(I)+1<=size(Air,I)
                xl(I)=xl(I)+1;
            end
        end
    end
end


