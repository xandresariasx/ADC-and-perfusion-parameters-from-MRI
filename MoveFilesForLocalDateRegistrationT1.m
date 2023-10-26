function MoveFilesForLocalDateRegistrationT1(Folder,Date,ImageTemplate)


mkdir([Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Date filesep])
mkdir([Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Date filesep 'Processed' filesep 'Registered' filesep])
copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep ImageTemplate],...
    [Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Date filesep ImageTemplate])
copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep ImageTemplate],...
    [Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Date filesep 'Processed' filesep 'Registered' filesep ImageTemplate])
J=1;
for Str={'T1','AUC','Ve','Vp','Kt'} 
    if exist([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep Str{1} filesep])~=0   % 9/15/2020
        copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep Str{1} filesep],...
            [Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Date filesep Str{1} filesep])
        copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep Str{1} filesep],...
            [Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Date filesep 'Processed' filesep 'Registered' filesep Str{1} filesep])
    end
    J=J+1;
end

copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep 'TumorSegmentationT1' filesep],...
    [Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Date filesep 'TumorSegmentationT1' filesep])
copyfile([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date filesep 'Processed' filesep 'Registered_Date' filesep 'TumorSegmentationT1' filesep],...
    [Folder filesep 'T1W' filesep 'DateRegisteredLocal' filesep Date filesep 'Processed' filesep 'Registered' filesep 'TumorSegmentationT1' filesep])
