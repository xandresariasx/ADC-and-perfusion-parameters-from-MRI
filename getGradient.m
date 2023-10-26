function [GRAD,DIV] = getGradient(m,h)
%function [GRAD,DIV] = getGradient(m,h)
%
%For a given image-size m a discrete gradient is set-up using
%short-forward differences and Neumann boundary-conditions.
%Image dimensions of 1D,2D and 3D are supported [1].
%
%For the resulting matrix GRAD and a given image I it holds
%
%       1D                  2D                      3D
%       GRAD*I(:) = I1      GRAD*I(:)=[I1;I2]       GRAD*I(:)=[I1;I2;I3] 
%
%where Ii is the derivative of I in dimension i.
%Note that if I is given on a cell-centered grid, GRAD*I(:) is given on a
%staggered grid
% 
%INPUT:
%    m - Image Dimension, usually m=size(I)
%    h - Voxel spacing
%    
%OUTPUT:
%    GRAD - Discrete gradient
%     DIV - Discrete divergence, DIV=-GRAD'
%
%
%REFERENCES:
% [1] Modersitzki, "FAIR: flexible algorithms for image registration". 
%     Vol. 6. SIAM, 2009.
%
%                                       (c)Constantin Sandmann, 23-Nov-2015
%                                                 http://mic.uni-luebeck.de
%                                     See LICENSE.txt for copyright issues
% ------------------------------------------------------------------------- 

    %if no input is given: show help
    if nargin == 0 
        clc
        help(mfilename) 
        return; 
    end

    %setup main variables
    D   = @(i) getDerivative1D(i,m,h);
    I   = @(i) speye(m(i));
    dim = nnz(m-1);
    
    
    %do the job
    switch dim
        case 1
            GRAD = D(1);
             DIV = -GRAD';
             
        case 2
            GRAD = [kron(I(2),D(1));
                    kron(D(2),I(1))];
             DIV = -GRAD';
             
        case 3
            GRAD = [kron(I(3),kron(I(2),D(1)));
                    kron(I(3),kron(D(2),I(1)));
                    kron(D(3),kron(I(2),I(1)))];
             DIV = -GRAD';        
        otherwise
            error('Dimensionality not supported');
    end




end




function D = getDerivative1D(i,m,h)
%1D derivative with Neumann BC

    %get setup some variables
    ni = m(i);    
    hi = h(i);
    
    %setup derivative matrix
    e          = ones(ni,1)*[-1,1];
    D          = spdiags(e,[0,1],ni,ni)./hi;
    D(end,end) = 0;
    
end
