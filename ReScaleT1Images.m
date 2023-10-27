function [VolsRS,Info,Cs]=ReScaleT1Images(FolderT1Vols, FolderMasks, RefTiss)

if ~exist('RefTiss','var')
    RefTiss=[];
end

Files=AdjustDirVariable(dir([FolderT1Vols 'T1W_Alpha=*']));
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

load([FolderMasks filesep 'Masks.mat'])

% Get Average Intensity per image for muscle and liver
for II=1:length(Labels)
    IsLable{II}=[]; 
    for I=1:length(VolsRT1)
        Vol=VolsRT1{I};        
        %IsLable{II}(I)=mean(Vol(logical(MaksPerLabel{II}) & Vol~=0)); 
        IsLable{II}(I)=median(Vol(logical(MaksPerLabel{II}) & Vol~=0));   %9/7/2020 use median instead of mean in case of unstable masks
    end
end

% Define signal Eq. (S) and SSE function f

S=@(Mo,t1,a,tr) (Mo.*((1-exp(-tr/t1))./(1-cos(a).*exp(-tr/t1))).*sin(a));

Sp=@(t1,a,tr) (((1-exp(-tr/t1))./(1-cos(a).*exp(-tr/t1))).*sin(a));

g=@(x)  10000*(tanh(x-15)+1);

f = @(Mo,C,Alphas,Trs,Is,T1)   sum((Is-S(C*Mo,T1,Alphas,Trs)).^2) +  sum(g(C))     ;


% Optimize scaling values

options = optimset('MaxFunEvals',20000,'MaxIter',20000);

TableT1L={'Tissue',                                  '1.5T',       '3T';...
           'tumor',                                     0.7,        0.8;...
           'Muscle',                                      1,        1.4;...
           'liver',                                    0.57,       0.8;...
           'kidney/RenalCortex',                        0.9,        1.1;...
           'fat',                                      0.34,       0.35;...
           'bone/BoneMarrow',                          0.54,        0.58;...
           'spleen',                                      1,        1.3;...
           'aorta/vein/blood/artery/BloodVessel',       1.4,        1.9};

       
 for II=1:length(Labels)
    Is=IsLable{II};
    C_M1{II}=[]; C_M2{II}=[];
    if ~any(isnan(Is))
        aux=find(cellfun(@(x) contains(Labels{II},x,'IgnoreCase',true) | ...
            contains(x,Labels{II},'IgnoreCase',true),TableT1L(2:end,1)));        
        if ~isempty(aux) %&& ~strcmpi(Labels{II},'tumor')
            if Info{1}{1}.MagneticFieldStrength<3
                T1L=TableT1L{aux+1,2}; 
            else
                T1L=TableT1L{aux+1,3};
            end
            [x,~,ExFlag] = fminsearch(@(x) f(x(1),x(2:end),AnglesD*pi/180,trS,Is,T1L),[2*max(Is), ones(1,numel(Is))],options);
            C_M1{II}=x(2:end);

            Mo=Is(1)/Sp(T1L,AnglesD(1)*pi/180,trS(1));
            C=Is(2:end)./(Mo*Sp(T1L,AnglesD(2:end)*pi/180,trS(2:end)));
            C=[1 C];
            C_M2{II}=C;
        end
    end
 end
 
% Get errors using scaling values 
ErrorsOri=[];
ft = fittype('(Mo*((1-exp(-tr/t1))/(1-cos(a)*exp(-tr/t1)))*sin(a))',...
            'dependent',{'y'},'independent',{'a','tr'});
for II=1:length(Labels)    
    ErrorsM1{II}=[];
    ErrorsM2{II}=[];   
    if ~any(isnan(IsLable{II}))
        [fitresultO, ~,~] = fit([(AnglesD*pi/180)', trS'], IsLable{II}', ft, 'StartPoint',[IsLable{II}(1), 1],'Lower',[0 0],'Upper', [Inf Inf],...
                                'Display','off');
        aux=find(cellfun(@(x) contains(Labels{II},x,'IgnoreCase',true) | ...
            contains(x,Labels{II},'IgnoreCase',true),TableT1L(2:end,1)));               
        if Info{1}{1}.MagneticFieldStrength<3           
            T1L=TableT1L{aux+1,2};             
        else
            T1L=TableT1L{aux+1,3};
        end
        ErrorsOri(II)=abs(fitresultO.t1-T1L)/T1L;
        for I=1:length(Labels)  
            ErrorsM1{II}(I)=Inf;
            ErrorsM2{II}(I)=Inf;
            if ~isempty(C_M1{I})
                [fitresultM1, ~,~] = fit([(AnglesD*pi/180)', trS'], (IsLable{II}./C_M1{I})', ft, 'StartPoint',[IsLable{II}(1)./C_M1{I}(1), 1],'Lower',[0 0],'Upper', [Inf Inf],...
                                    'Display','off');                            
                [fitresultM2, ~,~] = fit([(AnglesD*pi/180)', trS'], (IsLable{II}./C_M2{I})', ft, 'StartPoint',[IsLable{II}(1)./C_M2{I}(1), 1],'Lower',[0 0],'Upper', [Inf Inf],...
                                    'Display','off');                             
                ErrorsM1{II}(I)=abs(fitresultM1.t1-T1L)/T1L;
                ErrorsM2{II}(I)=abs(fitresultM2.t1-T1L)/T1L;
            end
        end
    end
end

% return scaling that gives minimum error
if ~isempty(RefTiss)
    aux=contains(Labels,RefTiss,'IgnoreCase',true);
    [val,L]=min([ErrorsOri(aux) ErrorsM1{aux}(aux) ErrorsM2{aux}(aux)]);
    Vols={VolsRT1, cellfun(@(x,y) x./y, VolsRT1, num2cell(C_M1{aux}),'UniformOutput',false), ...
        cellfun(@(x,y) x./y, VolsRT1, num2cell(C_M2{aux}),'UniformOutput',false)};
    VolsRS=Vols{L};
    Cs=[{ones(1,length(VolsRT1))} C_M1(aux) C_M2(aux)];
    Cs=Cs{L};
    Methods={'O','1', '2'};
    Method=Methods{L};
    ErrorOri=ErrorsOri(aux);
    save([FolderT1Vols 'ScalingData.mat'],'ErrorOri','val','Method','RefTiss')
    save([FolderT1Vols 'ScalingData2.mat'],'Cs','AnglesD','trS')
    disp(['Selected method: Method ' Method ' Ref. ' Labels{find(aux)} ...
        ' (' num2str(val*100,1) '%)']) 
    return
else
   [val,L]=min([mean(ErrorsOri) min(mean(cat(1,ErrorsM1{:}),1)) min(mean(cat(1,ErrorsM2{:}),1))]);
   aux={ErrorsOri, ErrorsM1, ErrorsM2};
   aux=aux{L};
   if L~=1
        [~,L2]=min(mean(cat(1,aux{:}),1));
   end
   Cs={{ones(1,length(VolsRT1))} C_M1 C_M2};
   Cs=Cs{L};
   if L~=1
        Cs=Cs{L2};
   end
   VolsRS=cellfun(@(x,y) x./y, VolsRT1, num2cell(Cs),'UniformOutput',false);
   Methods={'O','1', '2'};
   Method=Methods{L};
   ErrorOri=mean(ErrorsOri);
   if L~=1
       RefTiss=['Auto: ' Labels{L2}];   
   end
   save([FolderT1Vols 'ScalingData.mat'],'ErrorOri','val','Method','RefTiss')
   save([FolderT1Vols 'ScalingData2.mat'],'Cs','AnglesD','trS')
   disp(['Selected method: Method ' Method ' Ref. ' RefTiss ...
        ' (' num2str(val*100,2) '%)']) 
end
 
 
 
 
 
 
 
 
 
 

 

