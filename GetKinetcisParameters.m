function [Vp,Ve,Kt]=GetKinetcisParameters(WriteFolder,Folder,FolderMask)

aux=strsplit(Folder,filesep);
aux(cellfun(@isempty, aux))=[];
Patient=aux{end-4};
Date=aux{end-1};

Air=[];
if exist([strjoin(aux(1:end-3),filesep) filesep 'AssumeT1sC.mat'],'file')   % 4/6/2021 Use assumed T1
    load([strjoin(aux(1:end-3),filesep) filesep 'AssumeT1sC.mat'])
    if AssumeT1sC==1
        VolAvg=GetAvDCEImg(FolderMask);
        Air=VolAvg/max(VolAvg(:))<0.01;   % Air in the image
    end
end
if isempty(Air)
    load([Folder 'T1Data_NoCorrection.mat'])        % 2/8/2021 to not use AirTh
    Air=T1==0;
end

load([Folder 'Concentrations.mat'])
Cons=Cons_NW_NL;
Cond=Conds_NW_NL;

AIF=GetAIFManualAnnotations(WriteFolder,Folder,Cons,Cond,Times,Air);
AIF=[Times'; AIF];

AIF(2,:)=AIF(2,:)/(1-0.45);     %Hematrocrit correction, is it neccessary?

save([Folder filesep 'AIF.mat'], 'AIF')

% Shift AIF
figure, plot(AIF(1,:),AIF(2,:))
load([FolderMask filesep 'Masks.mat']);
L=find(contains(Labels,'Tumor','IgnoreCase',true),1);
[x,y,z]=ind2sub(size(MaksPerLabel{L}),find(MaksPerLabel{L}(:)));
ConsT=[];
for II=1:length(x)
    if all(squeeze([Cond(x(II),y(II),z(II),:)])) 
            ConsT=[ConsT; squeeze([Cons(x(II),y(II),z(II),:)])']; 
    end
end
hold on, plot(Times,nanmedian(ConsT))
[~,Center] = kmedoids(ConsT,1);
hold on, plot(Times,Center)

Shift=input('Enter # of steps to shift AIF: ');  
if Shift>0
    Times=Times(1+Shift:end);
    Times=Times-Times(1);
    AIFs = interp1(AIF(1,:),AIF(2,:),Times,'linear');
    clear AIF;
    AIF(1,:)=Times';
    AIF(2,:)=AIFs';
    Cons=Cons(:,:,:,1+Shift:end);
    Cond=Cond(:,:,:,1+Shift:end);
end
%%%%%

[x, y, z, n]=size(Cons);
n=n-1;                      % Number of time points
TimesSec=AIF(1,:);          % time points in s

if isempty(Cond)
    Conds=ones(x, y, z);
else
    Conds=sum(Cond,4)/size(Cond,4);
end
Conds=Conds & ~Air;
aux=Conds==1;
N=sum(aux(:)==1);           % Total number of voxels to calculate

IndP=find(Conds==1);        % Get indexes of voxels


% Get column 1 of matrix A given by the integral of AIF from 0 to every
% time point
At_1=[];
for I=2:size(AIF,2)
    At_1=[At_1 trapz(TimesSec(1:I),AIF(2,1:I))];
end

% Get column 3 of matrix A given by the AIF at every time point
At_3=AIF(2,2:end);   %-AIF(2,1);

% Get the columns 2 of matrix A given by - the integrals of the
% concentrations from 0 to every time point

aux=reshape(Cons,[x*y*z n+1]);
Cts=aux(IndP,:);
At_2=[];
for I=2:size(AIF,2)
    At_2(:,I-1)=trapz(TimesSec(1:I),Cts(:,1:I),2);
end
At_2=-At_2;

% Get C vectors given by the concentrations at every voxel
Cs=Cts(:,2:end);   %-Cts(:,1);

% for every voxel solve C\A to obtain b
% vp=b(3), ve=b(1)/b(2)-b(3), and Ktrans=b(1)-b(2)*b(3)
for I=1:length(IndP)     % change to parfor
    A=[];
    A=At_1';
    A(:,2)=At_2(I,:)';
    A(:,3)=At_3';
    C=Cs(I,:)';
%     mdl=fitlm(A,C,'y ~ x1 + x2 + x3 - 1','RobustOpts','on');
%     B=mdl.Coefficients{:,1};
    B=A\C;
    vp(I)=B(3);
    ve(I)=(B(1)/B(2))-B(3);
    Ktrans(I)=B(1)-B(2)*B(3);  
    
    % Fit quality measure
%     Cts_fit=[];
%     Error=[];
%     Cts_fit(1)=0;
%     for J=2:size(AIF,2)
%             Cts_fit(J)=Ktrans(I)*trapz(TimesSec(1:J),AIF(2,1:J).*exp(-B(2)*(TimesSec(J)-TimesSec(1:J))))+vp(I)*AIF(2,J);
%     end
%     Cts_fits(I,:)=Cts_fit;
%     Errors(I)=sqrt(sum((Cts_fit-Cts(I,:)).^2)/mean(Cts(I,:))^2);   
    %%%%   
end

% Build 3D matrices with the solutions
Vp=zeros(x, y, z);
Vp(IndP)=vp;

Ve=zeros(x, y, z);
Ve(IndP)=ve;

Kt=zeros(x, y, z);
Kt(IndP)=Ktrans;


figure('units','normalized','outerposition',[0 0 1 1]), colormap('gray')       
[S,~]=ReadDcmFolder4([Folder 'DCE_t=5' filesep]);
S=S{1};
SlicePos=floor(size(S,3)/2);
Maps={Kt,Vp,Ve};
Titles={'Kt','Vp','Ve'};
Ranges={[0 0.0050],[0 0.2], [0 1],[0 0.5]};
for I=1:length(Maps)
    subplot(1,3,I), 
    Overlay(S(:,:,SlicePos), real(Maps{I}(:,:,SlicePos)), 0.5, Ranges{I})    
    title(Titles{I},'Interpreter','none')
end
savefig([Folder filesep 'KineticMaps'])
saveas(gcf,[Folder filesep 'KineticMaps.jpeg']) 





function AIF=GetAIFManualAnnotations(WriteFolder,Folder,Cons,Conds,Times, Air)

aux=strsplit(Folder,filesep);
aux(isempty(aux))=[];
Patient=aux{find(contains(aux,'T1W'))-1};
Date=aux{find(contains(aux,'Registered'))+1};

% Conditions to discard some voxels
I1min=find(Times/60 < 1,1,'last');
Cons1min=Cons(:,:,:,2:I1min);   % before 1:I1min
Consdecay=Cons(:,:,:,I1min+1:end);  
Av1min=mean(Cons1min,4);
Avdecay=mean(Consdecay,4);
I2min=find(Times/60 > 2 &...
    Times/60 < 3 ,1,'first');
I3min=find(Times/60 > 2 &...
    Times/60 < 3 ,1,'last');
Cons2to3min=Cons(:,:,:,I2min:I3min);
NoiseCons=std(Cons2to3min,0,4);
VoxelsAifShape=Av1min>Avdecay & (Av1min-Avdecay)>NoiseCons;
aux=sum(Cons<0,4);
VoxelsPosCons=aux==0;
aux=sum(Cons>10,4);
VoxelsInrangeCons=aux==0;
Cond2=VoxelsAifShape & VoxelsPosCons & VoxelsInrangeCons;
%%%%%%%
Cond3=mean(Cons,4)>0.1;  % 10/19/2020
Cond4=std(Cons(:,:,:,find(Times>200,1):end),[],4)<0.05*max(Cons,[],4);   % 2/17/2021 before (0.05) it was a very strong condition
Cond2=Cond2 & Cond3 & Cond4 & ~Air;
%%%%%%

Concentrations=Cons;
Conditions=Conds;

load([WriteFolder Patient filesep 'T1W' filesep Date filesep 'Masks.mat']);

L=find(contains(Labels,'aorta','IgnoreCase',true) | contains(Labels,'blood','IgnoreCase',true)...
    | contains(Labels,'artery','IgnoreCase',true));
[x,y,z]=ind2sub(size(MaksPerLabel{L}),find(MaksPerLabel{L}(:)));
Cons=[];
for II=1:length(x)
    if all(squeeze([Conditions(x(II),y(II),z(II),:)])) && Cond2(x(II),y(II),z(II))
            Cons=[Cons; squeeze([Concentrations(x(II),y(II),z(II),:)])']; 
    end
end

if isempty(Cons)    % 2/17/2021  remove Cond4 in case no curve is selected
    Cond2=VoxelsAifShape & VoxelsPosCons & VoxelsInrangeCons &...
        ~Air & Cond3;
    Cons=[];
    for II=1:length(x)
        if all(squeeze([Conditions(x(II),y(II),z(II),:)])) && Cond2(x(II),y(II),z(II))
                Cons=[Cons; squeeze([Concentrations(x(II),y(II),z(II),:)])']; 
        end
    end
end


if isempty(Cons)
    Cons=[];
    for II=1:length(x)
        if all(squeeze([Conditions(x(II),y(II),z(II),:)])) && ~Air(x(II),y(II),z(II))
                Cons=[Cons; squeeze([Concentrations(x(II),y(II),z(II),:)])']; 
        end
    end
    AIF=median(Cons,1);
else    
    AIF=median(Cons,1);
end



function Vol=GetAvDCEImg(Folder)

Files=dir([Folder filesep 'Processed' filesep 'Registered' filesep]);
FileNames={Files(:).name};
aux=cellfun(@(x) contains(x,'DCE_t='), FileNames);
FileNames=FileNames(aux);
Vols=cellfun(@(x) ReadDcmFolder4([Folder filesep 'Processed' filesep 'Registered' filesep x filesep]),...
    FileNames);
Vols=cat(4,Vols{:});
Vol=mean(Vols,4);  






















