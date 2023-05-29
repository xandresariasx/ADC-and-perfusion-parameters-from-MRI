function MoveFilesForDateRegistrationDW(Folder,Date,LocalRegTemplate)

try,rmdir([Folder '\DWI\DateRegistered\' Date '\'],'s');end
mkdir([Folder '\DWI\DateRegistered\' Date '\'])
mkdir([Folder '\DWI\DateRegistered\' Date '\Processed\Registered\'])
copyfile([Folder '\DWI\Registered\' Date '\Local\' LocalRegTemplate],...
    [Folder '\DWI\DateRegistered\' Date '\' LocalRegTemplate])
copyfile([Folder '\DWI\Registered\' Date '\Local\' LocalRegTemplate],...
    [Folder '\DWI\DateRegistered\' Date '\Processed\Registered\' LocalRegTemplate])
[~,Infos]=ReadDcmFolder4([Folder '\DWI\Registered\' Date '\Local\' LocalRegTemplate '\']);
load([Folder '\DWI\Registered\' Date '\Local\ADC.mat'])
WriteDicomFolderV2(ADCav, Infos{1}, [Folder '\DWI\DateRegistered\' Date '\ADC\'],'ADC');   
copyfile([Folder '\DWI\DateRegistered\' Date '\ADC\'],...
    [Folder '\DWI\DateRegistered\' Date '\Processed\Registered\ADC\'])


MaskDW=load([Folder '\DWI\Registered\' Date '\Global\Masks.mat']);  
Labels=MaskDW.Labels;
[VolDW,InfosDW]=ReadDcmFolder4([Folder '\DWI\Registered\' Date '\Local\' LocalRegTemplate '\']);
VolDW=VolDW{1};
InfosDW=InfosDW{1};
FieldsTumorDW=find(  (contains(Labels, 'tumor', 'IgnoreCase',true)  |   cellfun(@(x) ~isempty(regexpi(x,'(t)+[0-9]')), Labels)) ...
                    & ~contains(Labels, 'core', 'IgnoreCase',true) & ~contains(Labels, 'viable', 'IgnoreCase',true));
MatMasksTumors=cat(4,MaskDW.MaksPerLabel{FieldsTumorDW});
MatMasksTumors=any(MatMasksTumors,4);
WriteDicomFolderV2(MatMasksTumors, InfosDW, [Folder '\DWI\DateRegistered\' Date '\TumorSegmentationDW\'],...
    'TumorSegmentation',1,0) 
copyfile([Folder '\DWI\DateRegistered\' Date '\TumorSegmentationDW\'],...
    [Folder '\DWI\DateRegistered\' Date '\Processed\Registered\TumorSegmentationDW\'])





