function [S,dS] = simulateSignalV3(T1,M0,flips,TR)
%function S = simulateSignal(T1,M0,flips,TR)
%
%For T1,M0,flips,TR MR Signals are simulated according to the signal
%equation
%
%   S = sin(alpha)*M0*( 1-exp(-TR/T1) )/( 1-cos(alpha)*exp(-TR/T1) )
%
%
%INPUT:
%    T1 - T1-relaxation times, size(T1)=[m1,m2,m3].
%    M0 - M0 magnetization, size(M0)=[m1,m2,m3].
% flips - flip angles in RAD, size(flips)=[k,1];
%    TR - repetition time
% 
%OUTPUT:
%     S - Simulated signal, size(S)=[m1,m2,m3,k].
%    dS - Derivative of the signal with respect to (T1,M0). 
%         For n := m1*m2*m3, size(dS) = [n*k,2*n].
%
%
%REFERENCES:
% [1] Zur, Stokar et al, "An analysis of fast imaging sequences with 
%     steady-state transverse magnetization refocusing," Magn Reson Med, 
%     6(2), 1988.
%
%                                      (c)Constantin Sandmann, 23-Nov-2015
%                                                http://mic.uni-luebeck.de
%                                     See LICENSE.txt for copyright issues
% ------------------------------------------------------------------------- 

    %if no input is given: show help
    if nargin == 0 
        clc
        help(mfilename) 
        return; 
    end
    
    
    %do Derivative only if the derivative is wanted
    doDerivative = (nargout>1);

    %prepare basic variables
    m     = size(T1);
    k     = size(flips,1);            %numel(flips);    
    n     = numel(T1);    
    %flips = reshape(flips,n,k)'; %reshape(flips,1,k);
    T1    = reshape(T1,n,1);
    M0    = reshape(M0,n,1);
    

    %setup helper variables
    E1  = exp(-TR./T1);
    c  = cos(flips);
    s  = sin(flips);
    ek = ones(1,k);
    en = ones(n,1);

    %simulate signal
%     S = (M0*s).*(1-E1*ek)./(1-E1*c);
    try       
        if size(E1,2)>1
           S = (repmat(M0,1,k).*s').*(1-E1)./(1-(E1.*c')); 
        else
            S = (repmat(M0,1,k).*s').*(1-E1)./(1-(repmat(E1,1,k).*c'));
        end
        S = squeeze(reshape(S,[m,k]));
    catch
        S = (M0*s).*(1-E1*ek)./(1-E1*c);
        S = reshape(S,[m,k]);
    end
    

    %if necessary: setup the derivative
    if doDerivative
        dE    = (TR./T1.^2).*E1;
        try
            if size(E1,2)>1
               dSdT  = (repmat(M0,1,k).*s').*(-dE)./(1-(E1.*c'))...
                    -(repmat(M0,1,k).*s').*(1-E1).*(-dE.*c')./...
                    (1-(E1.*c')).^2;
                dSdM  = (repmat(en,1,k).*s').*(1-E1)./(1-(E1.*c'));  
            else
                dSdT  = (repmat(M0,1,k).*s').*(-dE)./(1-(repmat(E1,1,k).*c'))...
                    -(repmat(M0,1,k).*s').*(1-E1).*(-repmat(dE,1,k).*c')./...
                    (1-(repmat(E1,1,k).*c')).^2;
                dSdM  = (repmat(en,1,k).*s').*(1-E1)./(1-(repmat(E1,1,k).*c'));   
            end
        catch
            dSdT  = (M0*s).*(-dE*ek)./(1-E1*c)...
                    -(M0*s).*(1-E1*ek).*(-dE*c)./(1-E1*c).^2;
            dSdM  = (en*s).*(1-E1*ek)./(1-E1*c);
        end        
        idxi  = (1:k*n);
        idxj  = repmat((1:n)',1,k);
        dS    = [sparse(idxi,idxj,dSdT,k*n,n), sparse(idxi,idxj,dSdM,k*n,n)];
    end

end