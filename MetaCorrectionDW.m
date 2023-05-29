function MetaCorrectionDW(Folder,Dates) 

for Date=Dates
    [~,InfADCo]=ReadDcmFolder4([Folder '\DWI\DateRegistered\' Date{1} '\ADC\']);
    [Vol,InfADCR]=ReadDcmFolder4([Folder '\DWI\DateRegistered\' Date{1} '\Processed\Registered_Date\ADC\']);
    for I=1:numel(InfADCR{1})
        InfADCR{1}{I}.RescaleIntercept=InfADCo{1}{1}.RescaleIntercept;
        InfADCR{1}{I}.RescaleSlope=InfADCo{1}{1}.RescaleSlope;  
        dicomwrite(uint32(Vol{1}(:,:,I)), InfADCR{1}{I}.Filename, InfADCR{1}{I},...
            'MultiframeSingleFile',false,'CreateMode','copy'); 
    end
end
