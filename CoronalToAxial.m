% load coronal image
[VolOrg,FilesOrg]=ReadDcmFolder2('\\hlm\data\home\4467692$\Desktop\536570\02_06_08\COR_STIR\');
VolOrg=VolOrg{1};

% rotate matrix to axial view
%VolOrgAxial=permute(VolOrg,[3 2 1]);
%VolOrgAxial=VolOrgAxial(end:-1:1,:,:);
VolOrgAxial=permute(VolOrg,[2 3 1]);
% VolOrgAxial=permute(VolOrgAxial,[2 1 3]);
% VolOrgAxial=VolOrgAxial(end:-1:1,end:-1:1,:);
% VolOrgAxial=permute(VolOrgAxial,[2 1 3]);


% change image properties
try
    rmdir('C:\Temp\Institution\Study\536570_2\02_06_08\COR_STIR_AX\','s')
end
mkdir('C:\Temp\Institution\Study\536570_2\02_06_08\COR_STIR_AX\')

zr=FilesOrg{1}{1}.SpacingBetweenSlices;
xryr=FilesOrg{1}{1}.PixelSpacing;
xr=xryr(1);
yr=xryr(2);

for I=1:size(VolOrgAxial,3)
    dicomwrite(uint16(VolOrgAxial(:,:,I)), ['C:\Temp\Institution\Study\536570_2\02_06_08\COR_STIR_AX\Im_'...
        num2str(I) '.dcm'],   'MultiframeSingleFile',false)
    
    info=dicominfo(['C:\Temp\Institution\Study\536570_2\02_06_08\COR_STIR_AX\Im_' num2str(I) '.dcm']);
    info.SpacingBetweenSlices=xr;
    info.SliceThickness=xr;
    info.PixelSpacing=[zr;yr];
    info.SliceThickness=xr;
    info.InstanceNumber=I;
    info.SeriesDescription='Cor Stir to Ax';
    dicomwrite(uint16(VolOrgAxial(:,:,I)), ['C:\Temp\Institution\Study\536570_2\02_06_08\COR_STIR_AX\Im_'...
        num2str(I) '.dcm'], info, 'MultiframeSingleFile',false,'CreateMode','copy')
end

nop=1;

