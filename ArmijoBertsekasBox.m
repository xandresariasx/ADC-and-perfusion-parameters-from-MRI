function [t,Yt,LSiter,LS] = ArmijoBertsekasBox(objFctn,Yc,dY,Jc,dJ,varargin)
%function [t,Yt,LSiter,LS] = ArmijoBertsekasBox(objFctn,Yc,dY,Jc,dJ,varargin)
%
%Implements the Armijo linesearch scheme [1] for box-constraints outlied 
%in [2]. See GaussNewtonProjected for more information.
%
%
%
%REFERENCES:
% [1] Nocedal, Wright, "Numerical optimization", Springer, 2006.
% [2] Bertsekas, "Projected Newton methods for optimization problems with 
%     simple constraints", SIAM J Control Optim, 20(2), 1982.
% [3] Modersitzki, "FAIR: flexible algorithms for image registration". 
%     Vol. 6. SIAM, 2009.
%
%
%Original function "Armijo" from FAIR by Jan Modersitzki [3].
%
%                          (c)Modified by Constantin Sandmann, 23-Nov-2015
%                                                http://mic.uni-luebeck.de
%                                     See LICENSE.txt for copyright issues
% ------------------------------------------------------------------------- 

    %help
    if nargin == 0, 
        clc;
        help(mfilename);
        return;
    end;

    LSMaxIter   = 10;           % max number of trials
    LSreduction = 1e-4;         % slope of line
    t           = 1;            % initial step
    proj        = @projection;

    for k=1:2:length(varargin), % overwrites default parameter
      eval([varargin{k},'=varargin{',int2str(k+1),'};']);
    end;

    if ~exist('idxR','var'); error('No idxR specified'); end;

    %compute the free indices
    idxF       = (1:length(Yc))';
    idxF(idxR) = [];

    L = [];

    for LSiter =1:LSMaxIter,
      Yt = Yc + t*dY;       % compute test value Yt
      Yt = proj(Yt);
      Jt = objFctn(Yt);      % evalute objective function


      FreeValue = t*dJ(idxF')*(dY(idxF));
      ResValue  = dJ(idxR')*(Yt(idxR)-Yc(idxR));
      if isempty(ResValue); ResValue = 0; end;


      LS = (  Jt <= Jc + LSreduction*( FreeValue + ResValue )  );

      if LS 
          break; 
      end;    % success, return
      L = [L;t,Jt];
      t = t/2;          % reduce t

    end;



    if LS, return; end;      % we are fine
    fprintf('Line Search failed - break\n');
    t = 0; Yt = Yc;        % take no action

end