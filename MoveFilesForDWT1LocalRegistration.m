function MoveFilesForDWT1LocalRegistration(Folder,Dates,LocalRegTemplateB,ImageTemplateB)

% Meta correction of registered images
for Date=Dates
    [~,InfADCo]=ReadDcmFolder4([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'ADC_' Date{1} filesep]);
    [Vol,InfADCR]=ReadDcmFolder4([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'Processed' filesep 'Registered_Date' filesep 'ADC_' Date{1} filesep]);
    for I=1:numel(InfADCR{1})
        InfADCR{1}{I}.RescaleIntercept=InfADCo{1}{1}.RescaleIntercept;
        InfADCR{1}{I}.RescaleSlope=InfADCo{1}{1}.RescaleSlope;  
        dicomwrite(uint32(Vol{1}(:,:,I)), InfADCR{1}{I}.Filename, InfADCR{1}{I},...
            'MultiframeSingleFile',false,'CreateMode','copy'); 
    end
end
% Move DW files
rmdir([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'Processed' filesep 'Registered' filesep],'s')
movefile([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'Processed' filesep 'Registered_Date' filesep],[Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'Processed' filesep 'Registered' filesep])
rmdir([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep LocalRegTemplateB],'s')
copyfile([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'Processed' filesep 'Registered' filesep LocalRegTemplateB],[Folder filesep 'T1-DW Registered' filesep 'DWI' filesep LocalRegTemplateB])
rmdir([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'TumorSegmentationDW' filesep],'s')
copyfile([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'Processed' filesep 'Registered' filesep 'TumorSegmentationDW' filesep],[Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'TumorSegmentationDW' filesep])
for Date=Dates
    rmdir([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'ADC_' Date{1} filesep],'s')
    copyfile([Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'Processed' filesep 'Registered' filesep 'ADC_' Date{1} filesep], [Folder filesep 'T1-DW Registered' filesep 'DWI' filesep 'ADC_' Date{1} filesep])
end
% Move T1 files
rmdir([Folder filesep 'T1-DW Registered' filesep 'T1W' filesep],'s')
copyfile([Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Dates{1} filesep ImageTemplateB],...
    [Folder filesep 'T1-DW Registered' filesep 'T1W' filesep ImageTemplateB])
copyfile([Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Dates{1} filesep ImageTemplateB],...
    [Folder filesep 'T1-DW Registered' filesep 'T1W' filesep 'Processed' filesep 'Registered' filesep ImageTemplateB])

% Generate T1 segmentation
% MaskT1=load([Folder filesep 'T1W' filesep Dates{1} filesep 'Masks.mat']);
% Labels=MaskT1.Labels;
% FieldsTumorT1=find(  (contains(Labels, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+[0-9]')), Labels)) ...
%                     & ~contains(Labels, 'core', 'IgnoreCase',true) & ~contains(Labels, 'viable', 'IgnoreCase',true));
% MatMasksTumors=cat(4,MaskT1.MaksPerLabel{FieldsTumorT1});
% MatMasksTumors=any(MatMasksTumors,4);
% [~,InfosT1]=ReadDcmFolder4([Folder filesep 'T1-DW Registered' filesep 'T1W' filesep ImageTemplateB filesep]);
% InfosT1=InfosT1{1};
aux=AdjustDirVariable(dir([Folder filesep 'T1W' filesep...
    'DateRegisteredLocal' filesep Dates{1} filesep 'Processed' filesep]));
copyfile([Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Dates{1} filesep 'Processed' filesep aux(end).name filesep 'TumorSegmentationT1' filesep],...
    [Folder filesep 'T1-DW Registered' filesep 'T1W' filesep 'TumorSegmentationT1' filesep])
% WriteDicomFolderV2(MatMasksTumors, InfosT1, [Folder filesep 'T1-DW Registered' filesep 'T1W' filesep 'TumorSegmentationT1' filesep],...
%     'TumorSegmentationT1',1,0)  
copyfile([Folder filesep 'T1-DW Registered' filesep 'T1W' filesep 'TumorSegmentationT1' filesep],...
    [Folder filesep 'T1-DW Registered' filesep 'T1W' filesep 'Processed' filesep 'Registered' filesep 'TumorSegmentationT1' filesep])





