function [r,Contrast]=Get_rContrast(Folder, Contrast)

[~,Info]=ReadDcmFolder4([Folder 'DCE_t=10' filesep]);
Info=Info{1};


try, MF=Info{1}.MagneticFieldStrength; catch, MF=1.5; end
if MF~=1.5 || MF~=3, MF=1.5; end



Table={'Agent',         'Field Strength',   'r';
       'Multihance',    1.5             ,   6.2;
       'Multihance',    3               ,  5.37;
       'Magnevist',     1.5             ,  4.25;
       'Magnevist',     3               ,  3.76;
       'Gadavist',      1.5             ,  4.61;
       'Gadavist',      3               ,  4.46;
       'Dotarem',       1.5             ,  3.91;
       'Dotarem',       3               ,  3.43;
       'Omniscan',      1.5             ,  4.47;
       'Omniscan',      3               ,  3.89;
       'Optimark',      1.5             ,  4.43;
       'Optimark',      3               ,  4.24;
       'Eovist',        1.5             ,  7.24;
       'Eovist',        3               ,  5.45;
       'ProHance',      1.5             ,  4.39;
       'ProHance',      3               ,  3.46};
   
try, 
    if contains(Info{1}.ContrastBolusAgent,Table(2:end,1),'IgnoreCase',true)  
        Contrast=Info{1}.ContrastBolusAgent; 
    end
end  
 for I=1:size(Table,1)
     if contains(Contrast,Table{I,1},'IgnoreCase',true)  
         if MF==Table{I,2}
             r=Table{I,3};
             break
         else
             if MF==Table{I+1,2}
                 r=Table{I+1,3};
                 break
             end
         end
     end
 end
 
fid=fopen([Folder filesep 'UsedContrast.txt'],'w');
try
fprintf(fid, [Contrast ' ' num2str(MF) ' ' num2str(r)]);
catch
    nop=1;
end
fclose(fid);
 
 
 