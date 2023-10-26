function [Cons,Conds, TimesSec]=GetConcentrations(Folder,FolderM, r,...
    ImageTemplate, ScaleMo, AssumeT1s, AssumeT1Artery)

if ~exist('ScaleMo','var')
    ScaleMo=0;
end

% Load Files
Files=AdjustDirVariable(dir([Folder 'DCE_t=*']));
ImagesNames={Files(:).name};
[~,Ind]=sort(cell2mat(cellfun(@(x) str2num(x(strfind(x,'=')+1:end)), ImagesNames,'UniformOutput' , false)));
Files=Files(Ind);
[Vols,Infs]=arrayfun(@(x) ReadDcmFolder4([x.folder filesep x.name filesep]), Files);
try    % 2/9/2021 because wrong metadata after registration
    cellfun(@(x) x{1}.RepetitionTime/1000, Infs);
catch
    CorrectTagsDicom(Folder);
    [Vols,Infs]=arrayfun(@(x) ReadDcmFolder4([x.folder filesep x.name filesep]), Files);
end

% Get C
Angles=GetAngles(Folder);
load([Folder filesep 'C30.mat'])
if numel(C30)>1
    C30=mean(C30(Infs{1}{1}.FlipAngle==Angles));
end

% Solve T1 equation
% syms tr t1 a Mo S positive
% assume(a>0 & a<pi/2)
% [t1, param, cond]=solve(S==Mo*((1-exp(-tr/t1))/(1-cos(a)*exp(-tr/t1)))*sin(a),t1, 'ReturnConditions', true);
t1 =@(S,Mo,a,tr) -tr./log((S - Mo*sin(a))./(S*cos(a) - Mo*sin(a)));

% Get T10 and Mo
if ~AssumeT1s
    Data=load([Folder filesep 'T1Data_NoCorrection.mat']); 
    Mo=Data.M0;  
    T1Map=Data.T1;   
    Mo=C30*Mo;   
else
    T1Map=AssumeT1([FolderM 'Masks.mat'],1.5);
    [Vol,Info]=ReadDcmFolder4([Folder 'DCE_t=1' filesep]);
    trS=Info{1}{1}.RepetitionTime/1000;
    AnglesD=Info{1}{1}.FlipAngle;
    AnglesR=AnglesD*pi/180; 
    A=((1-exp(-trS(1)/T1Map))./(1-cos(AnglesR(1)).*exp(-trS(1)/T1Map))).*sin(AnglesR(1));
    Mo=Vol{1}./ A;
end
if AssumeT1Artery  % 5/11/2021  In case T1s at artery are very off
    T1MapAux=AssumeT1([FolderM 'Masks.mat'],1.5);
    [Vol,Info]=ReadDcmFolder4([Folder 'DCE_t=1' filesep]);
    trS=Info{1}{1}.RepetitionTime/1000;
    AnglesD=Info{1}{1}.FlipAngle;
    AnglesR=AnglesD*pi/180; 
    A=((1-exp(-trS(1)/T1MapAux))./(1-cos(AnglesR(1)).*exp(-trS(1)/T1MapAux))).*sin(AnglesR(1));
    MoAux=Vol{1}./ A;
    load([FolderM 'Masks.mat'])
    aux=cell2mat(cellfun(@(x) find(contains(Labels,x,'IgnoreCase',true)),{'artery','blood','aorta'},...
        'UniformOutput',false));
    ArterySeg=MaksPerLabel{aux};
    Mo(ArterySeg==1)=MoAux(ArterySeg==1);
    T1Map(ArterySeg==1)=T1MapAux(ArterySeg==1);
end

% Solve T1 for each time point and compute concentration
load([FolderM 'Masks.mat'])
aux=cell2mat(cellfun(@(x) find(contains(Labels,x,'IgnoreCase',true)),{'artery','blood','aorta'},...
    'UniformOutput',false));
ArterySeg=MaksPerLabel{aux};
IsMed=cellfun(@(x) median(x(ArterySeg==1)), Vols);
CombineAcqusitionAndTriggerTimes=1;
% if  any(cellfun(@(x) contains(Folder,x),...
%             {'102-003-104','102-003-105','HZ001-001-301','HZ001-001-905','HZ001-003-306',...
%                 'HZ001-003-904','HZ001-007-901','HZ001-007-903','HZ001-007-904'}))
%     CombineAcqusitionAndTriggerTimes=0;
% end
if CombineAcqusitionAndTriggerTimes
    TrigTimes=cell2mat(cellfun(@(x) GetRealTime2(x{1}), Infs, 'UniformOutput', false));
else    
    TrigTimes=cellfun(@(x) x{1}.AcquisitionTime, Infs, 'UniformOutput', false);
    TrigTimes=cell2mat(cellfun(@(x) GetRealTime(x), TrigTimes, 'UniformOutput', false));    
end
if mean(TrigTimes(2:end)-TrigTimes(1:end-1))>100
    TrigTimes=TrigTimes/1000;
end
TimesSec=(TrigTimes-TrigTimes(1))';
Sigma=12;    % before 7   ,  9/10/2020
[Gder,GderN]=GaussDerivative(TimesSec,IsMed,Sigma);
%k=find(GderN>0.7,1);           % k=find(GderN>0.45,1);
[~,k]=max(GderN);

%%%%%
fig=figure; yyaxis left, plot(TimesSec,IsMed)
hold on, yyaxis left, plot(TimesSec(k),IsMed(k),'o')
hold on, yyaxis right, plot(TimesSec,GderN,'r')
%%%%

t = timer('ExecutionMode', 'singleShot', 'StartDelay', 600, 'TimerFcn', @pressEnter);
start(t);
res = input('Enter shift of injection point (o):');
if isempty(res)
    res=0;
end
stop(t);
delete(t);
k=k+res;
%%%%%
close(fig)
fig=figure; yyaxis left, plot(TimesSec,IsMed)
hold on, yyaxis left, plot(TimesSec(k),IsMed(k),'o')
hold on, yyaxis right, plot(TimesSec,GderN,'r')
%%%%

TrigTimes=TrigTimes(k:end);
TimesSec=TrigTimes-TrigTimes(1);

Ss=Vols(1:k);
trs=cellfun(@(x) x{1}.RepetitionTime/1000, Infs(1:k));
as=cellfun(@(x) x{1}.FlipAngle*pi/180, Infs(1:k));


tr=mode(trs);
a=mode(as);
Ss(trs~=tr | as~=a)=[];
Ss=cat(4,Ss{:});
S=mean(Ss,4);
So=S;
save([Folder filesep 'AverageImgBaseline.mat'], 'So');

if ScaleMo
    C=ReScaleMo(So, Mo, a, tr, T1Map, Folder);
    Mo=C*Mo;
end

if exist([Folder filesep 'AnglesR_corrected.mat'])
    load([Folder filesep 'AnglesR_corrected.mat'])
    a=AnglesR_corrected2{deg2rad(Angles)==a};
end

%%%  Instead of using Mo use T1map instead 5/20/2021
Mo=S./(((1-exp(-tr/T1Map))./(1-cos(a).*exp(-tr/T1Map))).*sin(a));
%%%

% T10av=eval(t1);
T10av=t1(S,Mo,a,tr);
T10=T10av;

C=[];
Conds=[];
C{1}=zeros(size(T10));
Conds{1}=ones(size(T10));
Conds{1}(imag(T10)>0 | isnan(T10))=0;
for I=k+1:length(Vols)     
    S=Vols{I};
    Info=Infs{I}{1};     
    tr=Info.RepetitionTime/1000;
    a=Info.FlipAngle*pi/180;
%     t1map=eval(t1);
    t1map=t1(S,Mo,a,tr);
    CondMap=~(imag(t1map)>0) &  ~isnan(t1map) & ~isinf(t1map); 
    C{end+1}=((1/t1map)-(1/T10))/r;
    Conds{end+1}=CondMap;    
end

Cons=cat(4,C{:});

if ~isempty(Conds)
    Conds=cat(4,Conds{:});
end
       
save([Folder filesep 'T10av.mat'],'T10av')


        
        
   
function C=ReScaleMo(So, Mo, a, tr, T1Map, Folder)


aux=dir(Folder);
aux2=aux(contains({aux(:).name},'T1W_Alpha='));
VolsI=arrayfun(@(x) ReadDcmFolder4([x.folder filesep x.name filesep]),aux2);
aux=cat(4,VolsI{:});
Conds=sum(aux~=0,4)==size(aux,4);
Conds(T1Map==0)=0;

T10AvSc=@(C,S,Mo,Alpha,Tr)  -Tr./log((S - C*Mo*sin(Alpha))./(S*cos(Alpha) - C*Mo*sin(Alpha))) ;
f= @(C,S,Mo,Alpha,Tr,Conds,T1Map)  sum(   (T10AvSc(C,S(Conds(:)),Mo(Conds(:)),Alpha,Tr)-  T1Map(Conds(:))).^2     )  ;

T10Av=T10AvSc(1,So,Mo,a,tr);
Conds(imag(T10Av)>0 | isnan(T10Av))=0;

options = optimset('MaxFunEvals',20000,'MaxIter',20000);
[C,~,ExFlag] = fminsearch(@(x) f(x, So, Mo, a, tr, Conds, T1Map),1,options);


function [Gder,GderN]=GaussDerivative(TimesSec,IsMed,Sigma)

DGauss=@(x,Sigma) (-x./(Sigma^3*sqrt(2*pi))).*exp(-(x.^2)/(2*(Sigma^2)));
DTime=TimesSec(2:end)-TimesSec(1:end-1);
DTime=[DTime DTime(end)];
for I=1:length(IsMed)
    accum=0;
    for J=1:length(IsMed)
        accum=accum+IsMed(J)*DGauss(TimesSec(I)-TimesSec(J),Sigma)*DTime(J);  % 5s sigma
    end
    Gder(I)=accum;
end

GderN=Gder/max(Gder);

            
function AnglesD=GetAngles(Folder)

Files=AdjustDirVariable(dir([Folder 'T1W_Alpha=*']));
VolsRT1=[];Info=[];
for I=1:length(Files)
    [VolsRT1{I},Info{I}]=ReadDcmFolder3([Folder Files(I).name filesep]);
    VolsRT1{I}=VolsRT1{I}{1}; Info{I}=Info{I}{1};
end
AnglesD=cellfun(@(x) x{1}.FlipAngle, Info);
[~,Order]=sort(AnglesD);
AnglesD=AnglesD(Order);
  


function T1=AssumeT1(FileMask,FS)

load(FileMask)
TableT1L={'Tissue',                      '1.5T',       '3T';...
           'tumor',                                     0.7,        0.8;...
           'Muscle',                                      1,        1.4;...
           'liver',                                    0.57,       0.8;...
           'kidney/RenalCortex',                        0.9,        1.1;...
           'fat',                                      0.34,       0.35;...
           'bone/BoneMarrow',                          0.54,        0.58;...
           'spleen',                                      1,        1.3;...
           'aorta/vein/blood/artery/BloodVessel',       1.4,        1.9};
       
T1=ones(size(MaksPerLabel{1}));
for I=1:numel(Labels)
    TI=find(cellfun(@(x) contains(Labels{I},x,'IgnoreCase',true) | ...
        contains(x,Labels{I},'IgnoreCase',true),TableT1L(2:end,1)))+1;
    FI=find(strcmp(TableT1L(1,2:3),[num2str(FS) 'T']),1)+1;
    T1v=TableT1L{TI,FI};
    T1(MaksPerLabel{I}==1)=T1v;
end
        
        



