function MetaCorrectionDW(Folder,Dates) 

for Date=Dates
    [~,InfADCo]=ReadDcmFolder4([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date{1} filesep 'ADC' filesep]);
    [Vol,InfADCR]=ReadDcmFolder4([Folder filesep 'DWI' filesep 'DateRegistered' filesep Date{1} filesep 'Processed' filesep 'Registered_Date' filesep 'ADC' filesep]);
    for I=1:numel(InfADCR{1})
        InfADCR{1}{I}.RescaleIntercept=InfADCo{1}{1}.RescaleIntercept;
        InfADCR{1}{I}.RescaleSlope=InfADCo{1}{1}.RescaleSlope;  
        dicomwrite(uint32(Vol{1}(:,:,I)), InfADCR{1}{I}.Filename, InfADCR{1}{I},...
            'MultiframeSingleFile',false,'CreateMode','copy'); 
    end
end
