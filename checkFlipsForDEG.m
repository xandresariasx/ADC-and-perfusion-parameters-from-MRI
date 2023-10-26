function checkFlipsForDEG(flips)
%function checkFlipsForDEG(flips)
%
%Checks if the flip angles are given in DEG or in RAD by simply comparing 
%if they are larger than pi or smaller than -pi. I have made this this 
%mistake too often.
%
%
%
%
%
%                                      (c)Constantin Sandmann, 23-Nov-2015 
%                                                http://mic.uni-luebeck.de
%                                     See LICENSE.txt for copyright issues
%------------------------------------------------------------------------- 

    if nnz(flips>pi) || nnz(flips<-pi)
        warning('Are your flip-angles really given in RAD? Please check your input.')
    end

end