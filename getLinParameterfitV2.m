function [T1,M0,para] = getLinParameterfitV2(data,flips,TR,varargin)
%function [T1,M0,para] = getLinParameterfit(data,flips,TR,varargin)
%
%For given data, sequence parameters and weights w \in \R^k a linear
%parameterfit is performed to recover T1 by solving the voxel-wise
%weighted optimization problem:
%
%   Minimize:
%     J(z):= 1/2*|A*z - b|_W^2
%          = 1/2*(Az-b)'*W*(Az-b)
%
%Here A,W,b are given by 
%
%     A := [ d1./tan(flip1) , 1 ] ,   b := [ d1./sin(flip1) ]
%          [    .             . ]          [        .       ]
%          [    .             . ]          [        .       ]
%          [ dk./tan(flipk) , 1 ]          [ dk./sin(flipk) ]       
%
%      W := [w1,0 ,..,0 ]
%           [0 ,w2,..,0 ]
%           [0 ,0, ..,wk]
%
%The problem is solved by using an explicit inverse of (A'*W*A) in the  
%correspondingof normal equation (A'*W*A)*x = (A'*W)*b.
%
%INPUT:
%    data - dataset with 2D or 3D images acquired with k different flip 
%           angles. It needs to hold :
%           size(data)=[m1,m2,k] or size(data)=[m1,m2,m3,k].
%   flips - flip-angles, size(flips)=[k,1]
%      TR - repetition time, T1 will have the same units as TR.
%
%DEFAULT PARAMETERS:
%        w - Weights to be used. Default are the ones proposed by Gupta[1].
%            DEFAULT: ones(k,1)
%    T1Min - lower bound for T1, 
%            DEFAULT: 1
%    T1Max - upper bound for T1
%            DEFAULT: 1e5
%    M0Min - lower bound for M0, 
%            DEFAULT: 0
%    M0Max - upper bound for M0
%            DEFAULT: Inf
%    For more input arguments see the default parameters in the code
%
%OUTPUT:
%     T1 - Recovered T1-Times, same unit as TR
%     M0 - Recovered M0-Times
%   para - Struct with various output parameters. Fields:
%           E1 - exp(-TR/T1);
%            N - M0*(1-E1);
%          res - The residue S(z)-d
%
%
%REFERENCES:
% [1] Gupta, "A new look at the method of variable nutation angle for the 
%     measurement of spin-lattice relaxation times using fourier transform 
%     NMR," J Magn Reson, 25(1), 1977.
% [2] Fram, Herfkens, Johnson, et al., "Rapid calculation of T1 using 
%     variable flip angle gradient refocused imaging," Magn Reson Med, 
%     5(3), 1987.
% [3] Sandmann, Hodneland, Modersitzki, "A Practical Guideline for T1 
%     Reconstruction from Various Flip Angles in MRI", submitted 11/2015
%
%                                      (c)Constantin Sandmann, 23-Nov-2015
%                                                http://mic.uni-luebeck.de
%                                     See LICENSE.txt for copyright issues
% ------------------------------------------------------------------------- 


    %if no input is given: Show help and run minimal example
    if nargin == 0
        clc;
        help(mfilename);
        E2_LinearReconstructionOpt;
        clearvars T1 M0 para;
        return;
    end
    
    %default parameters
    w     = ones(size(data,4),1);         %numel(flips),1);
    T1Min = 1;
    M0Min = 0;
    T1Max = 1e5;
    M0Max = Inf;

    %overwrites default parameter 
    for j=1:2:length(varargin),  
        eval([varargin{j},'=varargin{',int2str(j+1),'};']); 
    end 

    
    %prepare variables
    %checkFlipsForDEG(flips);
    mk    = size(data);
    m     = mk(1:end-1);
    k     = mk(end);
    n     = prod(m);
    if numel(size(flips))>2
        flips = reshape(flips,n,k)';
    else
        flips = reshape(flips,k,1);
    end
    data  = reshape(data,n,k)';
    w     = reshape(w,k,1);

    % Solve the normal equation for
    % 
    %    J(x):=|V*(A*x - b)|^2
    % 
    % for V:=sqrt(W) and A,b as explained above
    
    %setup matrix V
    v = repmat(sqrt(w),n,1);
    V = spdiags(v,0,k*n,k*n);
    
    %setup matrix A=[AE1,AN]
    e     = ones(k*n,1);
    datat = bsxfun(@rdivide,data,tan(flips));
   
    idxi = (1:k*n);
    idxj = repmat(1:n,k,1);
    AE1  = sparse(idxi,idxj,datat(:),k*n,n);
    AN   = sparse(idxi,idxj,e(:),k*n,n);
    A    = [AE1,AN];
    
    %setup vector b
    b = bsxfun(@rdivide,data,sin(flips));
    
    %solve the system and get z=[E1;N]
    lhs  = V*A;
    rhs  = V*b(:);
    z    = lhs\rhs;
    E1   = z(1:n);
    N    = z(n+1:end);
    
    %apply constraints
    TR=mode(TR);
    E1(E1<exp(-TR/T1Min)) = exp(-TR/T1Min);
    E1(E1>exp(-TR/T1Max)) = exp(-TR/T1Max);
  

    %get T1
    T1 = reshape(-TR./log(E1),m);
    M0 = reshape(N./(1-E1),m);
    
    %apply constraints
    M0(M0<M0Min) = M0Min;
    M0(M0>M0Max) = M0Max;
    

    %simulate current model
    S   = simulateSignalV3(T1,M0,flips,TR);
    data=data';                 % Andres Arias
    res = S(:) - data(:);
    E1  = exp(-TR./T1);
    N   = M0.*(1-E1);
    
    
    %some output
    para = struct('res',reshape(res,mk),'E1',reshape(E1,m),'N',reshape(N,m));
    
    

end