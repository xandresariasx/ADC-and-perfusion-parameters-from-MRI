function [T1,M0,para] = getNLParameterfitDiffV2(data,h,flips,TR,gammaT1,...
    gammaM0,Air,T10,M00,Woption,varargin)
%function [T1,M0,para] = getNLParameterfitDiff(data,h,flips,TR,gammaT1,gammaM0,varargin)
%
%For given data and sequence parameters flips and TR, a stabilized,
%constrained nonlinear parameterfit is performed to recover T1 by solving 
%the voxel-wise optimization problem [1]:
%
%   Minimize:
%
%     J(z):= 1/2*|S(z) - d|^2 
%            + gammaT1/2*|\nabla z1|^2 + \gammaM0/2*|\nabla z2|^2
%
%where S(z) is the signal equation [2] and \nabla is a discrete gradient
%operator.
%In this implementation results from the linear parameterfit are used as
%starting guess [3,4]. The optimization is performed using a constrained 
%Gauss-Newton Algorithm proposed in [5]. Default constraints for T1 and M0 
%are listed under DEFAULT PARAMETERS.
%
%INPUT:
%    data - dataset with 2D or 3D images acquired with k different flip 
%           angles. It needs to hold :
%           size(data)=[m1,m2,k] or size(data)=[m1,m2,m3,k].
%       h - voxel spacing
%   flips - flip-angles, size(flips)=[k,1]
%      TR - repetition time, T1 will have the same units as TR.
% gammaT1 - stabilization parameter for T1
% gammaM0 - stabilization parameter for M0
%
%DEFAULT PARAMETERS:
%    T1Min - lower bound for T1, 
%            DEFAULT: T1Min=1
%    T1Max - upper bound for T1
%            DEFAULT: T1Max=1e5
%    M0Min - lower bound for M0, 
%            DEFAULT: M0Min=0
%    M0Max - upper bound for M0
%            DEFAULT: M0Max=Inf
%    For more input arguments see the default parameters in the code
%
%OUTPUT:
%     T1 - Recovered T1-Times, same unit as TR
%     M0 - Recovered M0-Times
%   para - Struct with various output parameters. Fields:
%           E1 - exp(-TR/T1);
%            N - M0*(1-E1);
%          res - The residue (S(z) - d)
%
%
%REFERENCES:
% [1] de Pasquale, Sebastiani, et al., "Bayesian estimation of relaxation 
%     times T1 in MR images of irradiated Fricke-agarose gels," 
%     Magn Reson Imag, 18(6), 2000. 
% [2] Zur, Stokar et al, "An analysis of fast imaging sequences with 
%     steady-state transverse magnetization refocusing," Magn Reson Med, 
%     6(2), 1988.
% [3] Wang, Shi, et al. "Concatenated and parallel optimization for the 
%     estimation of T1 maps in FLASH MRI with multiple flip angles", 
%     Magn Reson Med, 63(5), 2010.
% [4] Gupta, "A new look at the method of variable nutation angle for the 
%     measurement of spin-lattice relaxation times using fourier transform 
%     NMR," J Magn Reson, 25(1), 1977.
% [5] Bertsekas, "Projected Newton methods for optimization problems with 
%     simple constraints", SIAM J Control Optim, 20(2), 1982.
%
%                                      (c)Constantin Sandmann, 23-Nov-2015
%                                                http://mic.uni-luebeck.de
%                                     See LICENSE.txt for copyright issues
% ------------------------------------------------------------------------- 

    %if no input is given: Show help and run minimal example    
    if nargin == 0
        clc;
        help(mfilename);
        E4_DiffusiveStabilization;
        clearvars T1 M0 para;
        return;
    end
    
    
    %default parameters
    tolJ      = 1e-6; %stopping criterium for obFun value
    tolY      = 1e-6; %stopping criterium for update
    tolG      = 1e-6; %stopping criterium for gradient
    tolCG     = 1e-4; %stopping criterium for CG-iteration
    maxIter   = 50;   %maximal number of iterations for Gauss-Newton
    maxIterCG = 1000; %maximal number of iterations for CG
    LSMaxIter = 30;   %maximal number of iterations for line-search
    solver    = 'MatlabInternal'; %solver for the optimization
    regH      = 0;    %regularization of Hessian
    T1Min     = 1;    %constraints for T1
    T1Max     = 5000;  %constraints for T1
    M0Min     = 0;    %constraints for M0    
    M0Max     = Inf;  %constraints for M0


    %overwrites default parameter 
    for j=1:2:length(varargin),  
        eval([varargin{j},'=varargin{',int2str(j+1),'};']); 
    end 
    
    %check if flips are given in DEG
    checkFlipsForDEG(flips);
    
    %setup initial guess
    if ~exist('T10','var')
        T10=[];
    end
    if isempty(T10)  %~exist('T10','var')
%         for I=1:size(data,4)
%             w(I)=1/(RicianSTD(data(:,:,:,I))^2);
%         end
        lfParam   = {'T1Min',T1Min,'T1Max',T1Max,'M0Min',M0Min,'M0Max',M0Max};%,'w',w};
        [T10,M00] = getLinParameterfitV2(data,flips,TR,lfParam{:});
%         [T10,M00] =NLT1mapping(data,flips,TR,Air,0);
%         figure, Overlay(data(:,:,floor(size(data,3)/2),1), T10(:,:,floor(size(data,3)/2)), 0.5, [0 3000]) 
%         title('T10')
    end
    delete(gcp('nocreate'));   % 1/8/21

    %prepare optimization
    param0   = [T10(:);M00(:)];    
    obFun    =@(param) obFunRegParameterfit(param,data,h,flips,TR,gammaT1,gammaM0,Air,Woption);
    optParam = {'tolJ',tolJ,'tolY',tolY,'tolG',tolG,'tolCG',tolCG,...
                'maxIter',maxIter,'LSMaxIter',LSMaxIter,'maxIterCG',maxIterCG,...
                'solver',solver,'regH',regH,...
                'T1Min',T1Min,'T1Max',T1Max,'M0Min',M0Min,'M0Max',M0Max};
            
    %do the optimization
    paramOpt = GaussNewtonProjected(obFun,param0,optParam{:});            
            
    %get results
    [~,paraObFun] = obFun(paramOpt);
    
    %main output
    T1 = paraObFun.T1;
    M0 = paraObFun.M0;
    
    %pipe out some variables
    para = struct('res',paraObFun.res,'E1',paraObFun.E1,'N',paraObFun.N);

end


function [J,para,dJ,H] = obFunRegParameterfit(param,data,h,flips,TR,gamma1,gamma2,Air,Woption)
%function [J,para,dJ,H] = obFunRegParameterfit(param,data,flips,TR,gamma1,gamma2)
%Calculates the L2-Distance and corresponding derivatives.
%
% J = |S(z)-d|^2 + gamma1*|\nabla z1|^2 + gamma1*|\nabla z1|^2
%

    %setup main variables
    doDerivative = (nargout>2);
    
    %prepare variables
    mk    = size(data);
    m     = mk(1:end-1);
    n     = prod(m);
    T1    = reshape(param(1:n),m);
    M0    = reshape(param(n+1:end),m);
    k     = mk(end);
    if numel(size(flips))>2
        flips = reshape(flips,n,k)';
    else
        flips = reshape(flips,k,1);
    end
    
    Air=zeros(size(T1));
    if doDerivative
        [D,paraD,dD,d2D] = getDistance(T1,M0,data,h,flips,TR,Air,Woption);
        [R,dR,d2R]       = getRegularizer(T1,M0,h,gamma1,gamma2,Air);
        
        J  =   D +   R;
        dJ =  dD +  dR;
        H  = d2D + d2R;        
    else
        [D,paraD] = getDistance(T1,M0,data,h,flips,TR,Air,Woption);
        R         = getRegularizer(T1,M0,h,gamma1,gamma2,Air);        
    
        J  =   D + R;
    end    
    
    
    %output
    E1   = []; %exp(-mean(TR)./T1);
    N    = []; %M0.*(1-E1);
    para = struct('res',reshape(paraD.res,prod(m),[]),'T1',reshape(paraD.T1,m),'M0',reshape(paraD.M0,m),'E1',[],'N',[]);


end



function [D,para,dD,d2D] = getDistance(T1,M0,data,h,flips,TR,Air,Woption)


    %doDerivative switch
    doDerivative = (nargout>2);
    
    %prepare variables
    mk    = size(data);
    m     = mk(1:end-1);
    hd    = prod(h);
    T1    = reshape(T1,m);
    M0    = reshape(M0,m);
    
    %get weights
    switch Woption
        case 1
            w=ones(1,size(data,4));  % No weighting
        case 2
            for I=1:size(data,4)    % Weighting inverse to normalized variance 
                w(I)=1/(RicianSTD(data(:,:,:,I))^2);
            end
            w=w./max(w);
        case 3          % Weighting = I/var(I)
           for I=1:size(data,4)
                aux=data(:,:,:,I);
                w(I)=sqrt(mean(aux(:)))/(RicianSTD(data(:,:,:,I)));        
            end
            w=w.^2;
    end
       
    w=arrayfun(@(x) x.*ones(prod(m),1),w,'UniformOutput',false);
    w=cat(1,w{:});
    
    %get distance and (if necessary) derivative
    if doDerivative
        
        %simulate signal
        [S,dS] = simulateSignalV3(T1,M0,flips,TR);
        
        %setup obFun etc.
        res  = S(:) - data(:);        
        dres = dS;
        
        D=hd/2*(w.*res)'*res;
        
        k=size(flips,1);
        n=prod(m);
        idxi  = (1:k*n);
        idxj  = repmat((1:n)',1,k);
%         MatW    = [sparse(idxi,idxj,sqrt(w),k*n,n), sparse(idxi,idxj,sqrt(w),k*n,n)];
        MatW    = [sparse(idxi,idxj,(w),k*n,n), sparse(idxi,idxj,(w),k*n,n)];  
        dD   = hd*(res'*(dres.*MatW));
%         d2D  = hd*(dres'*dres);
        d2D  = hd*(dres'*(dres.*MatW));


    else

        
        %simualate Signal
        S = simulateSignalV3(T1,M0,flips,TR);
        
        %setup obFun etc.
        res  = S(:) - data(:);                
        D=hd/2*(w.*res)'*res;
        %D    = hd/2*(res'*res);
    end
    
    para = struct('res',reshape(res,mk),'T1',reshape(T1,m),'M0',reshape(M0,m));
end



function [R,dR,d2R] = getRegularizer(T1,M0,h,gammaT1,gammaM0,Air)
    
    %doDerivative switch
    doDerivative = (nargout>1);
    m            = size(T1);
    hd           = prod(h);
    n            = prod(m);
    persistent gammaT1c gammaM0c StabT1 StabM0
    
    %setup regularization matrix if necessary
    if isempty(StabT1) || (gammaT1c~=gammaT1) || (size(StabT1,2)~=n)
        GRAD     = getGradient(m,h);
        Si       = GRAD'*GRAD;
        StabT1   = gammaT1*hd*Si;
        gammaT1c = gammaT1;
    end
    if isempty(StabM0) || (gammaM0c~=gammaM0) || (size(StabM0,2)~=n)
        GRAD     = getGradient(m,h);
        Si       = GRAD'*GRAD;
        StabM0   = gammaM0*hd*Si;
        gammaM0c = gammaM0;
    end    
    
    %prepare variables
    T1 = T1(:);
    M0 = M0(:);
    
    T1f=T1;%.*~Air(:);
    M0f=M0;%.*~Air(:);
    
    if doDerivative
        R   = 1/2*(T1f'*StabT1*T1f + M0f'*StabM0*M0f);
        dR  = [T1f'*StabT1,M0f'*StabM0];
        d2R = blkdiag(StabT1,StabM0);
    else        
        R   = 1/2*(T1f'*StabT1*T1f + M0f'*StabM0*M0f);
    end

end