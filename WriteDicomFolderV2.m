function WriteDicomFolderV2(Vol, Info, filename, Description,m,Offset)

mkdir(filename)

if nargin<5
    m=[];
    Offset=[];
end
if isempty(Description)
    aux=strsplit(filename,'\');
    Description=aux{end};
end
Vol=double(Vol);
% Vol(isnan(Vol))=min(Vol(:));   % 10/26/2020
% Vol(isinf(Vol))=min(Vol(:));
Vol(isnan(Vol))=min(Vol(~isinf(Vol(:))));
Vol(isinf(Vol))=min(Vol(~isinf(Vol(:))));
if isempty(Offset)
    Offset=-prctile(-Vol(Vol<0),80);
    if isnan(Offset)
        Offset=0;
    end
end
Offset=round(Offset,3);
Vol=Vol-Offset;
if isempty(m)    
    m=(2^32-1)/max(Vol(:));
    if m<1
        m=(2^32-1)/prctile(Vol(:),95);
    end
end

m=floor(m);
Vol=Vol*m;    

for J=1:size(Vol,3)    
    if iscell(Info)
        info=Info{J};
    else
        info=Info(J); 
    end
    info.SeriesDescription=Description; 
    info.RescaleIntercept=Offset;           
    info.RescaleSlope=1/m;   
    try
        dicomwrite(uint32(Vol(:,:,J)), [filename num2str(J) '.dcm'],...
            info, 'MultiframeSingleFile',false,'CreateMode','copy')
    catch
        try
            dicomwrite(uint32(Vol(:,:,J)), [filename num2str(J) '.dcm'],...
                info, 'MultiframeSingleFile',false,'CreateMode','Create')
        catch
            try
                dicomwrite(uint16(Vol(:,:,J)), [filename num2str(J) '.dcm'],...
                    info, 'MultiframeSingleFile',false,'CreateMode','copy')
            catch
                dicomwrite(uint16(Vol(:,:,J)), [filename num2str(J) '.dcm'],...
                    info, 'MultiframeSingleFile',false,'CreateMode','Create')
            end
        end
    end
end