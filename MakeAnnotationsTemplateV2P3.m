function MakeAnnotationsTemplateV2P3(Folder, I, Template)

clc
close all force
warning off 
global Stop
global Stop2
Colors={'r','g','y','m','c','r','g','y','m','c','r','g','y','m','c','r','g','y','m','c'};
try, delete('Positions.mat'); end


if exist('Template','var')~=0
    Vol=ReadDcmFolder4([Folder Template '\']);
    Vol=Vol{1};
    aux=strsplit(Template,'_');
    Template2=strjoin({'b=450' aux{2:end}},'_');
    try
        Vol2=ReadDcmFolder4([Folder Template2 '\']);
        Vol2=Vol2{1};
    end
else
    if exist([Folder '\Processed\Registered\'])~=0
        Files=dir([Folder '\Processed\Registered\']);
        FileNames={Files(:).name};
        aux=cellfun(@(x) contains(x,'DCE_t='), FileNames);
        if all(aux==0)  %9/14/2020
            aux=cellfun(@(x) contains(x,'T1W_Alpha=30_1'), FileNames);
        end
        FileNames=FileNames(aux);
        Vols=cellfun(@(x) ReadDcmFolder4([Folder '\Processed\Registered\' x '\']),...
            FileNames);
        Vols=cat(4,Vols{:});
        Vol=mean(Vols,4);    
    else
        Vol=ReadDcmFolder4([Folder 'DCEav\']);
        Vol=Vol{1};
    end 
    Template='Average';
    try
        aux=strsplit(Folder,'\');
        aux2=cellfun(@isempty, aux);
        aux(aux2)=[];
        load([strjoin(aux(1:end-2),'\') '\T1W\ImageTemplates.mat']);
        Vol2=ReadDcmFolder4([Folder '\Processed\Registered\' ImageTemplate{I} '\']);
        Vol2=Vol2{1};
        Template2=ImageTemplate{I};
    end
end

figure('Name',Template,'OuterPosition',[378   477   576   513]), imshow3DOverlayCont( Vol, [],[],Folder,floor(size(Vol,3)/2),[],[],[])
if exist('Vol2','var')
    figure('Name',Template2,'OuterPosition',[943   477   576   513]), imshow3DOverlayCont( Vol2, [],[],Folder,floor(size(Vol,3)/2),[],[],[])
end

aux=strsplit(Folder,'\');
aux2=cellfun(@isempty, aux);
aux(aux2)=[];
if ~contains(Folder,'T1W')
    disp(['DWI Patient: ' aux{end-4} ' Date:' aux{end-1}]);
else
    disp(['DCE Patient: ' aux{end-2} ' Date:' aux{end}]);
end

Res=[];
while ~strcmpi(Res,'S') && ~strcmpi(Res,'R') && ~strcmpi(Res,'E')
    Res=input('Skip (S), Redo (R) or Edit (E) segmentations ','s');
end

if strcmpi(Res,'S')
    return
end

if strcmpi(Res,'R')
    close all force
    Fig=figure('Name',Template,'OuterPosition',[378   477   576   513]);imshow3DV3(Vol)
    set (Fig, 'WindowKeyPressFcn', @KeyPressFcn);
    set (Fig, 'WindowButtonDownFcn', @ButtonPressFcn);
    axis equal
    if exist('Vol2','var')
        Fig2=figure('Name',Template2,'OuterPosition',[943   477   576   513]);imshow3DV3(Vol2)
        set (Fig2, 'WindowKeyPressFcn', @KeyPressFcn);
        set (Fig2, 'WindowButtonDownFcn', @ButtonPressFcn);
        axis equal
    end
    Labels=input('Enter Labels to Segment: ','s');
    Labels=strsplit(Labels);    
    MaksPerLabel=[];
    PositionsPerLabel=[];
    for K=1:length(Labels)
        Stop=0;
        disp(['Make annotations: ' Labels{K}])
        %bw=[];
        positions=[];
        position=zeros(20);
        BW=zeros(size(Vol));
        I=1;
        while Stop==0 %size(position,1)>10 || size(positions,1)<10
            figure(Fig)
            if exist('Vol2','var')
                figure(Fig2)
            end
            disp('Select Image');
            Stop2=0;
            while Stop2==0
                pause(1)
            end
            h=imfreehand;
            if isempty(h)
                break
            end
            position = getPosition(h);
            bw=createMask(h);
        %     position=wait(h);
            aux2=get(gcf,'children');
            aux3=arrayfun(@(x) strcmp(class(x),'matlab.ui.control.UIControl'),...
                aux2);
            aux2_2=aux2(aux3);
            aux3=find(arrayfun(@(x) strcmp(x.Style,'slider'), aux2_2));
            SlicePos=floor(aux2_2(aux3).Value);
            BW(:,:,SlicePos)=BW(:,:,SlicePos) | bw;
            position=[position ones(size(position,1),1)*SlicePos];
            positions{I} = position; 
            ColorK=Colors{K};
            save('Positions.mat','positions','ColorK')
            positions2=cellfun(@(x) GetPositions(x,SlicePos), positions,'UniformOutput',false);
        %     positions2=positions(positions(:,3)==SlicePos,:);
            for J=1:length(positions2)
                hold on, plot(positions2{J}(:,1),positions2{J}(:,2),'Color',Colors{K})
            end
            delete(h)
            if exist('Fig2','var')                
                if isequal(Fig,gcf)
                    figure(Fig2)                    
                else
                    figure(Fig)
                end
                for J=1:length(positions2)
                    hold on, plot([NaN,NaN],'Color',Colors{K})
                end
                if isequal(Fig,gcf)
                    figure(Fig2)                    
                else
                    figure(Fig)
                end
            end
            I=I+1;
        end
        MaksPerLabel{K}=BW;
        PositionsPerLabel{K}=positions;
        aux=get(gca,'children');
        aux2=arrayfun(@(x) contains(class(x),'Line'),...
                                aux);
        aux2=aux(aux2);                        
        for II=1:length(aux2)
            set(aux2(II),'XData',[])
            set(aux2(II),'YData',[])
        end
        if exist('Fig2','var') 
            if isequal(Fig,gcf)
                figure(Fig2)                    
            else
                figure(Fig)
            end
            aux=get(gca,'children');
            aux2=arrayfun(@(x) contains(class(x),'Line'),...
                                    aux);
            aux2=aux(aux2);                        
            for II=1:length(aux2)
                set(aux2(II),'XData',[])
                set(aux2(II),'YData',[])
            end
        end
        delete('Positions.mat')
    end
    close all
    save([Folder 'Masks.mat'], 'MaksPerLabel','PositionsPerLabel','Labels');
    figure('Name',Template,'OuterPosition',[378   477   576   513]), imshow3DOverlayCont( Vol, [],[],Folder,floor(size(Vol,3)/2),[],[],[])
    if exist('Fig2','var') 
        figure('Name',Template2,'OuterPosition',[943   477   576   513]), imshow3DOverlayCont( Vol2, [],[],Folder,floor(size(Vol,3)/2),[],[],[])
    end
    disp('Saved...Enter to Continue')
    pause
    return
end

if strcmpi(Res,'E')
    load([Folder 'Masks.mat'])
    disp('Labels: ')
    disp(Labels);
    Res2=[];
    while ~strcmpi(Res2,'Y') && ~strcmpi(Res2,'N') 
        Res2=input('Delete labels? (Y/N)','s');
    end
    if strcmpi(Res2,'Y')
        Labels2=input('Enter Labels to Delete: ','s');
        Labels2=strsplit(Labels2);
        aux=cellfun(@(x) strcmp(Labels,x), Labels2,'UniformOutput' , false);
        aux=sum(cat(1,aux{:}),1)>=1;
        Labels(aux)=[];
        PositionsPerLabel(aux)=[];
        MaksPerLabel(aux)=[];
        disp('Labels: ')
        disp(Labels);
    end
    Res2=[];
    while ~strcmpi(Res2,'Y') && ~strcmpi(Res2,'N') 
        Res2=input('Add labels? (Y/N)','s');
    end
    if strcmpi(Res2,'Y')
        Labels2=input('Enter New Labels to Segment: ','s');
        Labels2=strsplit(Labels2);
        Labels=[Labels Labels2];
    end
    close all
    Fig=figure('Name',Template,'OuterPosition',[378   477   576   513]);imshow3DV3(Vol)
    set (Fig, 'WindowKeyPressFcn', @KeyPressFcn);
    set (Fig, 'WindowButtonDownFcn', @ButtonPressFcn);
    axis equal   
    if exist('Vol2','var')
        Fig2=figure('Name',Template2,'OuterPosition',[943   477   576   513]);imshow3DV3(Vol2)
        set (Fig2, 'WindowKeyPressFcn', @KeyPressFcn);
        set (Fig2, 'WindowButtonDownFcn', @ButtonPressFcn);
        axis equal
    end
    for K=1:length(Labels)
        Stop=0;
        disp(['Make annotations: ' Labels{K}])
        %bw=[];
        if K<=length(PositionsPerLabel)
            positions=PositionsPerLabel{K};
            ColorK=Colors{K};
            save('Positions.mat','positions','ColorK')
%             CopyDWContour;     % Remove!!!
            aux2=get(Fig,'children');
            aux3=arrayfun(@(x) strcmp(class(x),'matlab.ui.control.UIControl'),...
                aux2);
            aux2_2=aux2(aux3);
            aux3=find(arrayfun(@(x) strcmp(x.Style,'slider'), aux2_2));
            SlicePos=floor(aux2_2(aux3).Value);
            positions2=cellfun(@(x) GetPositions(x,SlicePos), positions,'UniformOutput',false);
            for J=1:length(positions2)
                hold on, plot(positions2{J}(:,1),positions2{J}(:,2),'Color',Colors{K})
            end
            if exist('Vol2','var')
                aux2=get(Fig2,'children');
                aux3=arrayfun(@(x) strcmp(class(x),'matlab.ui.control.UIControl'),...
                    aux2);
                aux2_2=aux2(aux3);
                aux3=find(arrayfun(@(x) strcmp(x.Style,'slider'), aux2_2));
                SlicePos=floor(aux2_2(aux3).Value);
                positions2=cellfun(@(x) GetPositions(x,SlicePos), positions,'UniformOutput',false);
                for J=1:length(positions2)
                    hold on, plot(positions2{J}(:,1),positions2{J}(:,2),'Color',Colors{K})
                end
            end
            BW=MaksPerLabel{K};
        else
            positions=[];
%             CopyDWContour;     % Remove!!!
            BW=zeros(size(Vol));
        end
        position=zeros(20);        
        I=length(positions)+1;
        while Stop==0 
            figure(Fig)
            if exist('Vol2','var')
                figure(Fig2)
            end
            disp('Select Image');
            Stop2=0;
            while Stop2==0
                pause(1)
            end
            h=imfreehand;
            if isempty(h)
                break
            end
            position = getPosition(h);
            bw=createMask(h);
            aux2=get(gcf,'children');
            aux3=arrayfun(@(x) strcmp(class(x),'matlab.ui.control.UIControl'),...
                aux2);
            aux2_2=aux2(aux3);
            aux3=find(arrayfun(@(x) strcmp(x.Style,'slider'), aux2_2));
            SlicePos=floor(aux2_2(aux3).Value);
            BW(:,:,SlicePos)=BW(:,:,SlicePos) | bw;
            position=[position ones(size(position,1),1)*SlicePos];
            positions{I} = position; 
            ColorK=Colors{K};
            save('Positions.mat','positions','ColorK')
%             CopyDWContour;     % Remove!!!
            positions2=cellfun(@(x) GetPositions(x,SlicePos), positions,'UniformOutput',false);
            for J=1:length(positions2)
                hold on, plot(positions2{J}(:,1),positions2{J}(:,2),'Color',Colors{K})
            end
            delete(h)
            if exist('Fig2','var')                
                if isequal(Fig,gcf)
                    figure(Fig2)                    
                else
                    figure(Fig)
                end
                for J=1:length(positions2)
                    hold on, plot([NaN,NaN],'Color',Colors{K})
                end
                if isequal(Fig,gcf)
                    figure(Fig2)                    
                else
                    figure(Fig)
                end
            end
            I=I+1;
        end
        MaksPerLabel{K}=BW;
        PositionsPerLabel{K}=positions;
        aux=get(gca,'children');
        aux2=arrayfun(@(x) contains(class(x),'Line'),...
                                aux);
        aux2=aux(aux2);                        
        for II=1:length(aux2)
            set(aux2(II),'XData',[])
            set(aux2(II),'YData',[])
        end
        if exist('Fig2','var') 
            if isequal(Fig,gcf)
                figure(Fig2)                    
            else
                figure(Fig)
            end
            aux=get(gca,'children');
            aux2=arrayfun(@(x) contains(class(x),'Line'),...
                                    aux);
            aux2=aux(aux2);                        
            for II=1:length(aux2)
                set(aux2(II),'XData',[])
                set(aux2(II),'YData',[])
            end
        end
        delete('Positions.mat')
    end
    save([Folder 'Masks.mat'], 'MaksPerLabel','PositionsPerLabel','Labels');
    figure('Name',Template,'OuterPosition',[378   477   576   513]), imshow3DOverlayCont( Vol, [],[],Folder,floor(size(Vol,3)/2),[],[],[])
    if exist('Fig2','var') 
        figure('Name',Template2,'OuterPosition',[943   477   576   513]), imshow3DOverlayCont( Vol2, [],[],Folder,floor(size(Vol,3)/2),[],[],[])
    end
    disp('Saved...Enter to Continue')
    pause
    return  
end



function PositionsOut=GetPositions(PositionsIn,SlicePos)

PositionsOut=PositionsIn(PositionsIn(:,3)==SlicePos,:);



function KeyPressFcn (object, eventdata)

global Stop
keyPressed = eventdata.Key;
if strcmp(keyPressed,'return')
    aux=get(gca,'children');
    delete(aux(1))
    Stop=1;
end

function ButtonPressFcn (object, eventdata)

global Stop2
Stop2=1;



    



