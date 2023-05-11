function dir2=AdjustDirVariable(dir2)

I3=1;
while I3<=length(dir2)    
    if isequal(dir2(I3).name,'.') || isequal(dir2(I3).name,'..')
        dir2(I3)=[];  
    else
        I3=I3+1;
    end
end