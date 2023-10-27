
function T1mapCorrectFA(Folder,FolderMask,gammaT1, gammaM0,...
    solver,AirTh,Date,Dates,WriteFolder,PatientName,WeightT1)

% Load Data
T1map(Folder,FolderMask,[], 0, gammaT1, gammaM0,...
    0, AirTh, solver,Date,Dates,WriteFolder,PatientName,0,WeightT1)
load([Folder 'T1Data_NoCorrection.mat'])
Mask=load([FolderMask 'Masks.mat']);
[Vols,Info]=ReadVolsT1(Folder);

% Get true T1 values from masks at the middle of the image

TableT1L={'Tissue',                      '1.5T',       '3T';...
           'tumor',                         0.7,        0.8;...
           'Muscle',                          1,        1.4;...
           'liver',                        0.57,       0.8;...
           'kidney/RenalCortex',            0.9,        1.1;...
           'fat',                          0.34,       0.35;...
           'bone/BoneMarrow',              0.54,        0.58;...
           'spleen',                          1,        1.3;...
           'aorta/vein/blood/artery',       1.4,        1.9};       


% Make B1 mask and correct FA

TR=Info{1}{1}.RepetitionTime/1000;
trS=cellfun(@(x) x{1}.RepetitionTime, Info)/1000;   
A=zeros(size(Mask.MaksPerLabel{1}));

for I=1:numel(Mask.Labels)    
    aux=find(contains(TableT1L(2:end,1),Mask.Labels{I},'IgnoreCase',true));
    if ~isempty(aux) && ~strcmpi(Mask.Labels{I},'tumor')
        T1t=TableT1L{aux+1,2}; 
        A(logical(Mask.MaksPerLabel{I}))=...
            acos(exp(-TR/T1t)) ./  acos(exp(-TR./T1(logical(Mask.MaksPerLabel{I}))));
    end    
end

AnglesD=cellfun(@(x) x{1}.FlipAngle, Info);
AnglesR=AnglesD*pi/180; 

AnglesR_corrected=arrayfun(@(x) x*A, AnglesR,'UniformOutput',false);

AnglesR_corrected2=AnglesR_corrected;
for z=1:size(T1,3)
    aux=cellfun(@(x) x(:,:,z), AnglesR_corrected,'UniformOutput',false);
    flips=cell2mat(cellfun(@(x) median(x(x~=0)),...
        aux,'UniformOutput',false));
    if isnan(flips(1))
        flips=AnglesR;
    end
    for J=1:numel(AnglesR_corrected2)
        AnglesR_corrected2{J}(:,:,z)=flips(J);
    end    
end
for I=1:numel(AnglesR_corrected2)
    aux=smooth(AnglesR_corrected2{I}(1,1,:));
    for J=1:numel(aux)
        AnglesR_corrected2{I}(:,:,J)=aux(J);
    end
end

h=[Info{1}{1}.PixelSpacing' Info{1}{1}.SpacingBetweenSlices];
Air=zeros(size(T1));
%Air(A~=0 | logical(Mask.MaksPerLabel{1}))=0;
% gammaT1 = 100*5e-5;                     % 5e-5;
% gammaM0 = 5e-5; 
% solver    = 'MatlabInternal';
[T1,M0] = getNLParameterfitDiffV2(cat(4,Vols{:}),h,cat(4,AnglesR_corrected2{:}),...
    trS*1000,gammaT1,gammaM0,Air,[],[],WeightT1,'solver',solver);

T1=T1./1000;

VolAvg=GetAvDCEImg(FolderMask);
Air=VolAvg/max(VolAvg(:))<AirTh;   % Air in the image
Air=Air | ~all(cat(4,Vols{:}),4); 
T1(Air)=0;        % 9/22/2020
M0(Air)=0;

save([Folder filesep 'AnglesR_corrected.mat'],'AnglesR_corrected2')
save([Folder filesep 'T1Data_NoCorrection.mat'],'T1','M0')

C30=1;
save([Folder filesep 'C30.mat'],'C30')



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




         