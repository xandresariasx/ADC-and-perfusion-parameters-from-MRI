function MoveFilesForLocalDateRegistrationDW(Folder,Date,LocalRegTemplate)

try,rmdir([Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Date filesep],'s');end
mkdir([Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Date filesep])
mkdir([Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Date filesep 'Processed' filesep 'Registered' filesep])
copyfile([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep LocalRegTemplate],...
    [Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Date filesep LocalRegTemplate])
copyfile([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep LocalRegTemplate],...
    [Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Date filesep 'Processed' filesep 'Registered' filesep LocalRegTemplate])
copyfile([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep 'ADC' filesep],...
    [Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Date filesep 'ADC' filesep])
copyfile([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep 'ADC' filesep],...
    [Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Date filesep 'Processed' filesep 'Registered' filesep 'ADC' filesep])

copyfile([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep 'TumorSegmentationDW' filesep],...
    [Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Date filesep 'TumorSegmentationDW' filesep])
copyfile([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep 'TumorSegmentationDW' filesep],...
    [Folder filesep 'DWI' filesep 'DateRegisteredLocal' filesep Date filesep 'Processed' filesep 'Registered' filesep 'TumorSegmentationDW' filesep])





