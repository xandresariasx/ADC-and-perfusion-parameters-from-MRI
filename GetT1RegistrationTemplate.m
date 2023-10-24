function ImageTemplate=GetT1RegistrationTemplate(FolderAux)

close all force
CombineAcqusitionAndTriggerTimes=1;

Images=AdjustDirVariable(dir([FolderAux 'DCE_t=*']));
if isempty(Images)
    ImageTemplate='T1W_Alpha=30_1';
    return
end
ImagesNames={Images(:).name};
[~,Ind]=sort(cell2mat(cellfun(@(x) str2num(x(strfind(x,'=')+1:end)), ImagesNames,'UniformOutput' , false)));
Images=Images(Ind);
[Vols,Infs]=arrayfun(@(x) ReadDcmFolder4([x.folder '\' x.name '\']), Images);

if CombineAcqusitionAndTriggerTimes
    TrigTimes=cell2mat(cellfun(@(x) GetRealTime2(x{1}), Infs, 'UniformOutput', false));
else    
    TrigTimes=cellfun(@(x) x{1}.AcquisitionTime, Infs, 'UniformOutput', false);
    TrigTimes=cell2mat(cellfun(@(x) GetRealTime(x), TrigTimes, 'UniformOutput', false));    
end


[Ut,~,It]=unique(TrigTimes);
for I=1:max(It)
    if sum(I==It)>1
        aux=cell2mat(cellfun(@(x) GetRealTime(x{1}.AcquisitionTime), Infs(I==It),...
            'UniformOutput', false)); 
        TrigTimes(I==It)=Ut(I)+aux-max(aux);
    end    
end

load([FolderAux 'Masks.mat']);
aux=cell2mat(cellfun(@(x) find(contains(Labels,x,'IgnoreCase',true)),{'artery','blood','aorta'},...
    'UniformOutput',false));
ArterySeg=MaksPerLabel{aux};
aux=cellfun(@(x,y) isequal(size(x),size(ArterySeg)), Vols);
Vols2=Vols(aux);

IsMed=cellfun(@(x) median(x(ArterySeg==1)), Vols2);
[~,k]=max(IsMed);

TimesSec=(TrigTimes-TrigTimes(1))';
TimesSec2=TimesSec(aux);


Titles=arrayfun(@num2str, [1:length(Vols2)],'UniformOutput' , false);
VolsHalf=cellfun(@(x) x(:,:,floor(size(x,3)/2)), Vols2,'UniformOutput' , false);
% MakeMultipleFigures(VolsHalf,Titles,[],[]);

figure('Name','Intensity Vs. time'), 
plot(TimesSec2,IsMed,'*-'), hold on, plot(TimesSec2(k),IsMed(k),'*r') 
saveas(gcf,[FolderAux '\MeanIntensityVsTime.jpeg'])

k=find(cumsum(aux)==k);

ImageTemplate=Images(k+2).name;

clc
close all force


