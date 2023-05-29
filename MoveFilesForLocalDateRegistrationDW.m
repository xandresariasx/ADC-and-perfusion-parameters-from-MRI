function MoveFilesForLocalDateRegistrationDW(Folder,Date,LocalRegTemplate)

try,rmdir([Folder '\DWI\DateRegisteredLocal\' Date '\'],'s');end
mkdir([Folder '\DWI\DateRegisteredLocal\' Date '\'])
mkdir([Folder '\DWI\DateRegisteredLocal\' Date '\Processed\Registered\'])
copyfile([Folder '\DWI\DateRegistered\' Date '\Processed\Registered_Date\' LocalRegTemplate],...
    [Folder '\DWI\DateRegisteredLocal\' Date '\' LocalRegTemplate])
copyfile([Folder '\DWI\DateRegistered\' Date '\Processed\Registered_Date\' LocalRegTemplate],...
    [Folder '\DWI\DateRegisteredLocal\' Date '\Processed\Registered\' LocalRegTemplate])
copyfile([Folder '\DWI\DateRegistered\' Date '\Processed\Registered_Date\ADC\'],...
    [Folder '\DWI\DateRegisteredLocal\' Date '\ADC\'])
copyfile([Folder '\DWI\DateRegistered\' Date '\Processed\Registered_Date\ADC\'],...
    [Folder '\DWI\DateRegisteredLocal\' Date '\Processed\Registered\ADC\'])

copyfile([Folder '\DWI\DateRegistered\' Date '\Processed\Registered_Date\TumorSegmentationDW\'],...
    [Folder '\DWI\DateRegisteredLocal\' Date '\TumorSegmentationDW\'])
copyfile([Folder '\DWI\DateRegistered\' Date '\Processed\Registered_Date\TumorSegmentationDW\'],...
    [Folder '\DWI\DateRegisteredLocal\' Date '\Processed\Registered\TumorSegmentationDW\'])





