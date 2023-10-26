function VOI=GenerateVOIforLocalDateRegistration(Folder,MaskFolder, ImageFolder, ImageTemplate, roiFactor,Infos)

load([MaskFolder filesep 'Masks.mat'])
FieldsTumor=find(  (contains(Labels, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+[0-9]')), Labels)) ...
                    & ~contains(Labels, 'core', 'IgnoreCase',true) & ~contains(Labels, 'viable', 'IgnoreCase',true));  
RGBs={[1 0 0],[0 1 0],[0 1 1],[1 1 0]};
cont=1;
VOI=[];
if ~isempty(ImageFolder)
[~,Infos]=ReadDcmFolder4([ImageFolder filesep 'Processed' filesep 'Registered' filesep ImageTemplate filesep]);
Infos=Infos{1}{1};
end
Resolution=[Infos.PixelSpacing' Infos.SliceThickness];
MatLabels=cat(4,MaksPerLabel{:});
aux2=[];
for I=1:length(MaksPerLabel),aux2{I}=I*ones(size(MaksPerLabel{1})); end
aux2=cat(4,aux2{:});
MatLabels=MatLabels.*aux2;
Dimensions=[size(MatLabels,1) size(MatLabels,2) size(MatLabels,3)];
if numel(roiFactor)==1
    roiFactor=ones(size(FieldsTumor))*roiFactor;
end
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
            Siz=Siz+max([0 0 0; 4-(Loc+Siz)]);   % 6/29/2020  To make voi greater than 4 in all dimensions
            Pos_mm=(Dimensions.*Resolution)./(Dimensions-1).*Loc-Dimensions.*Resolution.*...
                (0.5+(1./(Dimensions-1)));  
            Pos_mm=Pos_mm+Resolution/2;
            Siz_mm=Siz.*Resolution;
            roi=[Pos_mm(2) Pos_mm(1) Pos_mm(3) Siz_mm(2) Siz_mm(1) Siz_mm(3)];
            VOI(cont).xyzDxDyDz_mm=roi;
            VOI(cont).label=[Labels{FieldsTumor(I)} '_' num2str(cont)];
            VOI(cont).RGB=RGBs{mod(cont-1,4)+1};
            VOI(cont).index=cont;
            cont=cont+1;
        end
    end
end


