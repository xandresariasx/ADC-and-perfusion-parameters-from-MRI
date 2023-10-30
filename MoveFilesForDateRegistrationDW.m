function MoveFilesForDateRegistrationDW(Folder,Date,LocalRegTemplate)

try,rmdir([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep],'s');end
mkdir([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep])
mkdir([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep])
copyfile([Folder filesep 'DWI' filesep 'Registered' filesep Date filesep 'Local' filesep LocalRegTemplate],...
    [Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep LocalRegTemplate])
copyfile([Folder filesep 'DWI' filesep 'Registered' filesep Date filesep 'Local' filesep LocalRegTemplate],...
    [Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep LocalRegTemplate])
[~,Infos]=ReadDcmFolder4([Folder filesep 'DWI' filesep 'Registered' filesep Date filesep 'Local' filesep LocalRegTemplate filesep]);
load([Folder filesep 'DWI' filesep 'Registered' filesep Date filesep 'Local' filesep 'ADC.mat'])
WriteDicomFolderV2(ADCav, Infos{1}, [Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'ADC' filesep],'ADC');   
copyfile([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'ADC' filesep],...
    [Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep 'ADC' filesep])


MaskDW=load([Folder filesep 'DWI' filesep 'Registered' filesep Date filesep 'Global' filesep 'Masks.mat']);  
Labels=MaskDW.Labels;
[VolDW,InfosDW]=ReadDcmFolder4([Folder filesep 'DWI' filesep 'Registered' filesep Date filesep 'Local' filesep LocalRegTemplate filesep]);
VolDW=VolDW{1};
InfosDW=InfosDW{1};
FieldsTumorDW=find(  (contains(Labels, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+[0-9]')), Labels)) ...
                    & ~contains(Labels, 'core', 'IgnoreCase',true) & ~contains(Labels, 'viable', 'IgnoreCase',true));
MatMasksTumors=cat(4,MaskDW.MaksPerLabel{FieldsTumorDW});
MatMasksTumors=any(MatMasksTumors,4);
WriteDicomFolderV2(MatMasksTumors, InfosDW, [Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'TumorSegmentationDW' filesep],...
    'TumorSegmentation',1,0) 
copyfile([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'TumorSegmentationDW' filesep],...
    [Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered' filesep 'TumorSegmentationDW' filesep])





