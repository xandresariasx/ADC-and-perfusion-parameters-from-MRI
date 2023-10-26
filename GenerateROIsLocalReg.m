function obj=GenerateROIsLocalReg(obj,Infos, Labels, MatLabels, IncludeArtery, roiFactor) 

if exist('IncludeArtery','var')==0
    IncludeArtery=0;
end

% if roiFactor~=0     % 1/28/2021 to allow no local registration at tumor
FieldsTumor=find(  (contains(Labels, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+[0-9]')), Labels)) ...
                    & ~contains(Labels, 'core', 'IgnoreCase',true) & ~contains(Labels, 'viable', 'IgnoreCase',true));  
% else
%     FieldsTumor=[];
% end
if numel(roiFactor)==1
    roiFactor=ones(size(FieldsTumor))*roiFactor;
end
if IncludeArtery(1)
    FieldsArtery=find(contains(Labels,'artery','IgnoreCase',true) | contains(Labels,'aorta','IgnoreCase',true) | ...
        contains(Labels,'blood','IgnoreCase',true));
    FieldsTumor=[FieldsTumor FieldsArtery];
    if numel(IncludeArtery)>1
        roiFactor=[roiFactor IncludeArtery(2)];
    else
        roiFactor=[roiFactor 1.5];
    end
end

RGBs={[1 0 0],[0 1 0],[0 1 1],[1 1 0]};
cont=1;
obj.StudyDate.VOI_mm(:)=[];
obj.StudyDate.VOI(:)=[];
Resolution=[Infos.PixelSpacing' Infos.SliceThickness];
Dimensions=[size(MatLabels,1) size(MatLabels,2) size(MatLabels,3)];
for I=1:length(FieldsTumor)
    aux=MatLabels(:,:,:,FieldsTumor(I));
    if roiFactor(I)~=0
        CC=bwconncomp(aux);    
        for II=1:length(CC.PixelIdxList)
            aux2=zeros(size(aux));
            aux2(CC.PixelIdxList{II})=1;
            Stats=regionprops(aux2);

            Loc=round(Stats.BoundingBox(1:3));
            Siz=Stats.BoundingBox(4:6)-1;
            Loc=[Loc(2) Loc(1) Loc(3)];
            Siz=[Siz(2) Siz(1) Siz(3)];

            Factor=4.5./(Siz+1);
            Factor(roiFactor(I)*(Siz+1)-1>=3.5)=roiFactor(I);
            Loc=(Loc+Siz/2)-(Siz+1).*Factor/2;
            Loc=ceil(Loc); 
            Siz=Factor.*(Siz+1)-1;   

            [Loc,Siz]=FixVOI(obj,Loc,Siz);  

            if ~isempty(Siz)   
                Siz=Siz+max([0 0 0; 4-(Loc+Siz)]);        
                Siz=max([3 3 3; Siz]);  
                Loc=Loc+min([0 0 0; size(aux)-(Loc+Siz)]);

                Pos_mm=(Dimensions.*Resolution)./(Dimensions-1).*Loc-Dimensions.*Resolution.*...
                    (0.5+(1./(Dimensions-1)));  
                Pos_mm=Pos_mm+Resolution/2;

                Siz_mm=Siz.*Resolution;

                roi=[Pos_mm(2) Pos_mm(1) Pos_mm(3) Siz_mm(2) Siz_mm(1) Siz_mm(3)];


                obj.StudyDate.VOI_mm(cont,:)=roi;
                obj.StudyDate.VOI(cont).xyzDxDyDz_mm=roi;
                obj.StudyDate.VOI(cont).label=[Labels{FieldsTumor(I)} '_' num2str(cont)];
                obj.StudyDate.VOI(cont).RGB=RGBs{mod(cont-1,4)+1};
                obj.StudyDate.VOI(cont).index=cont;
                cont=cont+1;
            end
        end
    end
end




function [LocN,SizN]=FixVOI(obj,Loc,Siz)

aux={obj.StudyDate.vsRegistered(:).data};
aux=cat(4,aux{:});
aux2=all(aux,4);
for III=1:size(aux2,3)
    aux4=aux2(:,:,III);
    if sum(~aux4(:))>0.9*numel(aux4)
        aux2(:,:,III)=0;
    else
        aux2(:,:,III)=1;
    end
end
aux3=zeros(size(aux2));
aux3(max(Loc(1),1):min(Loc(1)+Siz(1),size(aux2,1)),...
    max(Loc(2),1):min(Loc(2)+Siz(2),size(aux2,2)),...
    max(Loc(3),1):min(Loc(3)+Siz(3),size(aux2,3)))=1;
aux4=and(aux2,aux3);
ind=find(aux4,1);
[x,y,z]=ind2sub(size(aux4),ind);
ind=find(aux4,1,'last');
[x2,y2,z2]=ind2sub(size(aux4),ind);
LocN=[x,y,z];
SizN=[x2-x,y2-y,z2-z];












