function AUC_W=Get90sAUC(Folder,FolderMask)

load([Folder 'Concentrations.mat'])
Cons=Cons_NW_NL;
Cond=Conds_NW_NL;

TimesSec=Times;

if isempty(Cond)
    CondT=ones(size(Cons,1),size(Cons,2),size(Cons,3));
else
    CondT=floor(sum(Cond,4)/size(Cond,4));
end
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

CondT=CondT & ~Air;    

pxi=linspace(0,TimesSec(end)-TimesSec(1),floor((TimesSec(end)-TimesSec(1))/0.5))';
pxi=pxi(pxi<90);
%pi=interp1((TimesSec-TimesSec(1)),permute(Cons,[4 1 2 3]),pxi,'linear');
%11/16/2020
[xi,xl]=GetCroppingInd(Air);
SizeO=size(Air);
clear T1 M0 Cons_NW_NL Conds_NW_NL Cond Air
Cons=Cons(xi(2):xl(2),xi(1):xl(1),xi(3):xl(3),:);
aux=permute(Cons,[4 1 2 3]); 
aux=real(aux);
%pi2=zeros(size(pxi,1),size(Cons,1),size(Cons,2),size(aux,4));
clear Cons 
for I=1:size(aux,4)
    aux2=aux(:,:,:,I);
    %aux2=gpuArray(aux2);     %  remove comment!
    pi=interp1((TimesSec-TimesSec(1)),aux2,pxi,'linear');
    %pi2(:,:,:,I)=gather(pi);   %  remove comment!
    pi2(:,:,:,I)=pi;    % Delete
    clear pi
end
pi=pi2;
clear pi2 aux2 aux
AUC_Wc=squeeze(trapz(pxi,pi,1));
AUC_W=zeros(SizeO);
AUC_W(xi(2):xl(2),xi(1):xl(1),xi(3):xl(3))=AUC_Wc;
AUC_W(~CondT)=0;



function Vol=GetAvDCEImg(Folder)

Files=dir([Folder filesep 'Processed' filesep 'Registered' filesep]);
FileNames={Files(:).name};
aux=cellfun(@(x) contains(x,'DCE_t='), FileNames);
FileNames=FileNames(aux);
Vols=cellfun(@(x) ReadDcmFolder4([Folder filesep 'Processed' filesep 'Registered' filesep x filesep]),...
    FileNames);
Vols=cat(4,Vols{:});
Vol=mean(Vols,4);  

function [xi,xl]=GetCroppingInd(Air)
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
[xi,xl]=CorrectXiXl(xi,xl,Air); 


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
















