function ADCMAPP3V2(Folder)

disp('Processing ADC all directions')
Files=AdjustDirVariable(dir([Folder '*']));
Files=Files([Files(:).isdir]);

if length(Files)<1
    return;
end

aux=cellfun(@(x) strsplit(x,'_'), {Files(:).name}, 'UniformOutput', false);
ScanNumFiles=cellfun(@(x) [x{3} '_' x{4}], aux,'UniformOutput',false);
ScanNum=unique(ScanNumFiles);
II=1;
for Scan=ScanNum
    Scan=Scan{1};
    FilesScanNum=Files(strcmp(ScanNumFiles,Scan));
    
    aux=cellfun(@(x) getb(x), {FilesScanNum(:).name}, 'UniformOutput', false);
    FilesScanNum(cellfun(@isempty,aux))=[];

    x=cell2mat(cellfun(@(x) getb(x), {FilesScanNum(:).name}, 'UniformOutput', false))';

    if length(unique(x))<=1    
        ADC{II}=[];
    else
        Vol=[];Info=[];    % 11/30/2020  this wasn't initialized before
        for I=1:length(FilesScanNum)
            [Vol{I},Info{I}]=ReadDcmFolder3([Folder FilesScanNum(I).name '\']);
            Vol{I}=Vol{I}{1}; Info{I}=Info{I}{1};
        end
        Vols=cat(4,Vol{:});
        ADC{II}=zeros(size(Vol{1}));
        So=zeros(size(Vol{1}));        
        for K=1:size(Vol{1},3)
            for J=1:size(Vol{1},2)
                for I=1:size(Vol{1},1)            
                    y=squeeze(Vols(I,J,K,:));
                    if all(y~=0)        
                        ly=log(y);
                        bs=[ones(length(ly),1) x]\ly;
                        ADC{II}(I,J,K)=-bs(2);
                        So(I,J,K)=exp(bs(1));
                    end
                end
            end
        end        
    end
    II=II+1;
end
ADC(cellfun(@isempty,ADC))=[];
ADCav=cat(4,ADC{:});
ADCav=mean(ADCav,4);
save([Folder '\ADC.mat\'],'ADC','ADCav')




 





function str_out=getb(str_in)


try
    aux=strsplit(str_in,'=');
    aux=strsplit(aux{2},'_');
    str_out=str2num(aux{1});
catch
    str_out=[];
end









