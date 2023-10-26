function SlicePos=GetSlicePos(Labels,MaksPerLabel,N)

MatLabels=cat(4,MaksPerLabel{:});
aux2=[];
for I=1:length(MaksPerLabel),aux2{I}=I*ones(N); end
aux2=cat(4,aux2{:});
MatLabels=MatLabels.*aux2;

FieldsTumor=find(contains(Labels,'tumor','IgnoreCase',true));
if isempty(FieldsTumor)            
    FieldsTumor=find(cellfun(@(x) ~isempty(regexpi(x,'(t)+[0-9]')), Labels));
end  
for I=1:N(3)
    aux=squeeze(MatLabels(:,:,I,:));
    aux5(I)=numel(unique(aux(:)))-1;
    aux6(I)=any(arrayfun(@(x) any(aux(:)==x), FieldsTumor));    
end
aux5(aux6==0)=0;
if sum(aux5(3:end-2)~=0)>0
    aux5(1:2)=0; aux5(end-1:end)=0;
end
SlicePos=find(aux5==max(aux5),1,'last');