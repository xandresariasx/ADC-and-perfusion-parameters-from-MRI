function [Vols,Files]=ReadDcmFolder3(Folder)

Dir=dir(Folder); 

Dir=AdjustDirVariable(Dir);

for I=1:length(Dir)
    aux{I}=Dir(I).name;
end

IDs=[];
for I=1:length(aux)
    aux2{I}=dicominfo([Folder aux{I}]);    
    IDs{end+1}= aux2{I}.SeriesInstanceUID;    
end

IDs=unique(IDs);

% get volumes ids and separate volumes
aux3=cell(length(IDs),1);
for I=1:length(IDs)
    for J=1:length(aux2)
        if isequal(aux2{J}.SeriesInstanceUID,IDs{I})
            aux3{I}{end+1}=aux2{J};
            %aux2{J}.InstanceNumber
        end  
    end
end

% sort images per Instance Number
for I=1:length(aux3)
    IN=cellfun(@(x) x.InstanceNumber, aux3{I});
    [~,IdIN]=sort(IN);
    aux4{I}=aux3{I}(IdIN);
end
aux3=aux4;

% Now create images from separated files
Im=cell(length(IDs),1);
for I=1:length(aux3)     
    for J=1:length(aux3{I})
        Aux=dicomread(aux3{I}{J});
        Im{I}{J}=squeeze(Aux);  
    end
end
for I=1:length(Im) 
    Im2=[];
    try
         for J=1:length(Im{I})
            Im2(:,:,J)=Im{I}{J};         
         end   
    catch
        Im2=Im{I}{1};
    end
     Vols{I}=Im2;
end

Files=aux3;

% It seems to work like this without checking the instance number


%imshow3D(Vols{I})





% 
% 
% for I=1:length(aux)
%     Aux=dicomread([Folder aux{I}]);
%     aux2=dicominfo([Folder aux{I}]);
%     try
%         if aux2.InstanceNumber==0 && length(aux)==1
%             Im=squeeze(Aux); 
%             return;
%         else 
%             if aux2.InstanceNumber==0
%                 Im=[];
%                 return;
%             end
%         end
%         if aux2.InstanceNumber~=I
%             nop=1;
%             Im=[];
%             return;
%         end
%     catch
%         Im=[];
%         return;
%     end
%     
%     Im{I}=squeeze(Aux);   
%                                     
%     %Im(:,:,I)= Aux;         
% end
% 
% aux=[];
% for I=1:length(Im)
%     aux=[aux ;size(Im{I})];    
% end
% 
% if length(unique(aux(:,1)))==1 && length(unique(aux(:,2)))==1 && size(aux,2)==2
%     for I=1:length(Im)
%         Im2(:,:,I)=Im{I};
%     end
%     Im=Im2;
% end









