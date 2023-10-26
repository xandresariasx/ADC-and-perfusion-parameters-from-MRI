function MetaCorrectionT1(Folder,Dates) 

for Date=Dates
    for Str={'T1','AUC','Ve','Vp','Kt'}
        try  % 9/15/2020
            [~,InfADCo]=ReadDcmFolder4([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date{1} filesep Str{1} filesep]);
            [Vol,InfADCR]=ReadDcmFolder4([Folder filesep 'T1W' filesep 'DateRegistered' filesep Date{1} filesep 'Processed' filesep 'Registered_Date' filesep Str{1} filesep]);
            for I=1:numel(InfADCR{1})
                InfADCR{1}{I}.RescaleIntercept=InfADCo{1}{1}.RescaleIntercept;
                InfADCR{1}{I}.RescaleSlope=InfADCo{1}{1}.RescaleSlope;  
                dicomwrite(uint32(Vol{1}(:,:,I)), InfADCR{1}{I}.Filename, InfADCR{1}{I},...
                    'MultiframeSingleFile',false,'CreateMode','copy'); 
            end
        end
    end
end
