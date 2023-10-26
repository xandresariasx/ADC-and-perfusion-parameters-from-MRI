function xc = GaussNewtonProjected(fctn,x0,varargin)
%function xc = GaussNewtonProjected(fctn,x0,varargin)
% 
%Gauss-Newton [1] optimization for T1 reconstruction with hard constraints 
%according to [2]. For x=[T1;M0] the problem
%
%    Minimize
%            J(x) := fctn(x)
%             s.t. T1Min<=T1<=T1Max and M0Min<=M0<=M0Max
%
%is solved.
% 
%INPUT:
%    fcnt - function handle for the objective function. 
%           For xc = [T1;M0] and size(xc)=[2*n,1] it must hold
%           [J,para,dJ,H] = fctn(x0), where:
%                J - is the function value, size(J)=[1,1]
%             para - struct with additional output parametes, may be empty.
%               dJ -  Jacobian, size(dJ) = [2*n,1]
%                H -  approximation to the Hessian, size(H) = [2*n,2*n].
%                     Hessian. Might be a function handle (this option is
%                     not testet, though).
%      x0 - Starting guess, x0 = [T10;M00] and size(x0)=[2*n,1].
%
%DEFAULT PARAMETERS:
%   T1Min - Minimal Value for T1
%           DEFAULT: 1
%   T1Max - Maximal Value for T1
%           DEFAULT: 1e4
%   M0Min - Minimal Value for M0
%           DEFAULT: 0
%   M0Max - Maximal Value for M0
%           DEFAULT: Inf
%   For more options see the default parameters in the code
%
%OUTPUT:
%      xc - Result of the optimization
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
%Original function "GaussNewton" from FAIR by Jan Modersitzki [3].
%
%                          (c)Modified by Constantin Sandmann, 23-Nov-2015
%                                                http://mic.uni-luebeck.de
%                                     See LICENSE.txt for copyright issues
% ------------------------------------------------------------------------- 

%if no input is given: Show help and run minimal example
if nargin ==0
  clc;
  help(mfilename);  
  runMinimalExample;
  return;
end;

% default parameters
maxIter      = 15;              % maximum number of iterations
maxIterCG    = 500;             % maximum number of iterations for CG
maxIterLS    = 30;              % maximum number of line search iterations
tolJ         = 1e-3;            % for stopping, objective function
tolY         = 1e-2;            %   - " -     , current value
tolG         = 1e-2;            %   - " -     , norm of gradient
tolCG        = 1e-4;            % stopping criterium for CG
solver       = [];              % solver to use
T1Min        = 1;               % constraints
T1Max        = 1e4;
M0Min        = 0;
M0Max        = Inf;
epsBdry      = 1e-5;            % distance to the boundary
idxF         = (1:numel(x0));   % free indices
idxR         = [];              % constrained indices
regH         = 0;               % Levenberg-Marquard parameter

for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%setup the projection function and vector-valued constraints
n       = numel(x0)/2;
constrL = [T1Min*ones(n,1);M0Min*ones(n,1)];
constrH = [T1Max*ones(n,1);M0Max*ones(n,1)];
proj    = @(x)  projection(x,constrL,constrH);



%initialize variables
xc                = proj(x0);
[JStop,para,dJ,H] = fctn(xc);

xOld = xc;
Jc   = JStop;
JOld = Jc;
iter = 0;
STOP = false(5,1);

%OUTPUT
fprintf('--------------------------------------------------------\n');
fprintf([mfilename, ' (CMS/JM 11/10/2015)\n']);
fprintf('--------------------------------------------------------\n\n');

fprintf('[ maxIter=%i / tolJ=%1.1e / tolG=%1.1e / tolY=%1.1e / length(yc)=%i ]\n\n',...
  maxIter,tolJ,tolG,tolY,numel(x0));
fprintf('iter\tJ\t\tJOld-Jc\t\t|\\nabla J|\tdy\t\tLS\n')
fprintf('%i\t%1.4e\t%1.4e\t%1.4e\t%1.4e\t%i\n',0,Jc,0,norm(dJ),0,0);       


%-- start the iteration --------------------------------------------------
while 1, 
   
  
    %get the diagonal of H to obtain epsc
    if isnumeric(H)
      di = spdiags(H,0);
    else
      try di = H.diag;
      catch, di = ones(numel(xc),1);
      end    
    end
  
    %compute updated epsilon_k
    wk   = norm(xc - proj(xc-(dJ')./di));
    epsc = min(epsBdry,wk);


    %get restricted indices
    idxL1 = find(xc <= constrL + epsc);
    idxL2 = find(dJ' > 0);
    idxRL = intersect(idxL1,idxL2);
    idxH1 = find(xc>=constrH - epsc);
    idxH2 = find(dJ' < 0);
    idxRH = intersect(idxH1,idxH2);
    idxR = union(idxRL,idxRH);
  
  
  
    %modify hessian  
    if isnumeric(H)
      %matrix which has only the diagonal of H
      Hdi     = spdiags(di,0,numel(xc),numel(xc));
      Htilde  = H;

      if ~isempty(idxR);
        Htilde(idxR,:) = Hdi(idxR,:);
        Htilde(:,idxR) = Hdi(:,idxR);
      end

      if regH ~= 0, Htilde = Htilde+regH*speye(length(xc));end
    else
      Htilde.fun  = @(x) multHtilde(x,H.fun,idxR,H.diag);      
      Htilde.diag = H.diag;
    end


    % solve the Gauss-Newton System
    dx = solveGN(-dJ',Htilde,solver,maxIterCG,tolCG);

  
    % check descent direction
    % note: descent is not granted if using an iterative solver 
    descent =   dJ * dx; 
    if descent > 0,
    warning('no descent direction, switch to -dy!')
    dx = -dx;
    end;

    % perform line-search
    [t,yt,LSiter] = ArmijoBertsekasBox(fctn,xc,dx,Jc,dJ,'para',para,...
        'LSMaxIter',maxIterLS,'idxR',idxR,'proj',proj);
    if (t == 0),
      warning('LineSearch failed.');
      break;
    end;


    % save old values and update
    xOld = xc; JOld = Jc; xc = yt; idxF = setdiff((1:length(xc)),idxR);
    iter = iter+1;

    [Jc,para,dJ,H] = fctn(xc); % evalute objective function

    %output
    fprintf('%i\t%1.4e\t%1.4e\t%1.4e\t%1.4e\t%i\n',iter,Jc,JOld-Jc,norm(dJ),norm(dx),LSiter);         

    %check stopping rules
    STOP(1) = (iter>0) && abs(JOld-Jc)   <= tolJ*(1+abs(JStop));
    STOP(2) = (iter>0) && (norm(xc-xOld) <= tolY*(1+norm(xc)));
    STOP(3) = norm(dJ(idxF))             <= tolG*(1+abs(JStop));
    STOP(4) = norm(dJ(idxF))             <= 1e6*eps;
    STOP(5) = (iter >= maxIter);
    if all(STOP(1:3)) || any(STOP(4:5)), break;  end;
  
end


%final output, if optimization is done
fprintf('STOPPING:\n');
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(1),...
  '(Jold-Jc)',abs(JOld-Jc),'tolJ*(1+|Jstop|)',tolJ*(1+abs(JStop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(2),...
  '|yc-xOld|',norm(xc-xOld),'tolY*(1+norm(yc)) ',tolY*(1+norm(xc)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(3),...
  '|dJ|',norm(dJ(idxF)),'tolG*(1+abs(Jstop))',tolG*(1+abs(JStop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(4),...
  'norm(dJ)',norm(dJ(idxF)),'eps',1e6*eps);
fprintf('%d[ %-10s=  %-14d >= %-25s=  %-14d]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);
fprintf('Number of restricted variables: %i of %i (ca. %2.2f %%) \n',length(idxR),length(xc),length(idxR)/length(xc)*100)
fprintf('-------------------done!--------------------------------\n');

end

function dx = solveGN(rhs,H,solver,maxIterCG,tolCG)
   
if isempty(solver) 
    if isnumeric(H),
      dx = H\rhs;
    else
      [dx,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG);
    end
else
    switch solver
        % matrix based
        % ------------
        case 'pcg'
            L   = tril(H); % Symmetric Gauss Seidel Preconditioning,
            D   = diag(H); % L is lower, D is diagonal, U = L'
            SGS = @(x) L\(D.*(L'\x));
            [dx,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,SGS);
        case 'ichol'
            %choose alpha such that H + alpha*diag(diag(H)) is diagonally dominant
            alpha = max(sum(abs(H),2)./diag(H))-2; 
            %setup options for ichol
            opts.diagcomp = alpha;
            opts.type = 'nofill';
            L1   = ichol(sparse(H),opts); % incomplete cholesky
            [dx,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,L1,L1');
        case 'jacobi-pcg'
            D   = diag(H); % D is diagonal
            PC = @(x) D.\x;
            [dx,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,PC);
        case 'cg'
%             [dx,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG);            
            try
                gpuDevice(1);
                Hgpu=gpuArray(H);
                [dx,flag,relres,iter] = pcg(Hgpu,rhs,tolCG,maxIterCG,Hgpu);
                dx=gather(dx);
            catch
                [dx,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,H);
            end 
        case 'MatlabInternal'
            dx = H\rhs;
        % matrix free
        case {'jacobi-pcg-mf'}
            D   = H.diag; % D is diagonal
            fun = H.fun;
            PC = @(x) x./D;
            [dx,flag,relres,iter] = pcg(fun,rhs,tolCG,maxIterCG,PC);
        case {'cg-mf'}
            fun = H.fun;
            [dx,flag,relres,iter] = pcg(fun,rhs,tolCG,maxIterCG);
        otherwise
            error('Solver not supported')
    end
end

if exist('flag','var')    
    switch flag
        case 1
            fprintf('pcg iterated %d times but converged only to relres %e instead of %e\n',...
                iter,relres,tolCG);
        case 2
            fprintf('Preconditioner M was ill-conditioned.\n');
        case 3
            fprintf('pcg stagnated. (Two consecutive iterates were the same.)\n');
        case 4
            fprintf('One of the scalar quantities calculated during pcg became too small or too large to continue computing.\n');
        otherwise
            fprintf('pcg success! %d iterations / relres= %1.2e / tolCG= %1.2e\n',iter,relres,tolCG);                    
    end
end




end


function x = projection(x,constrL,constrH)
%function xProj = projection(x,constrL,constrH)

    %do the projection
    x(x<constrL) = constrL(x<constrL);
    x(x>constrH) = constrH(x>constrH);

end




function Hx = multHtilde(x,Hfun,idx,diag)
    if isempty(idx)
        Hx = Hfun(x);
    else
        %save restricted elements and diagonal
        tmp = x(idx);
        d   = diag(idx);

        %H*(decoupled system)
        x(idx) = 0;
        Hx = Hfun(x);

        %rewrite indices
        Hx(idx) = tmp.*d;
    end
end

function runMinimalExample

    
    fprintf('====================MinimalExample======================\n');
    
    %setup an example functional J(x) := (x-b)'*(x-b) for size(x)=[n,1]
    n = 10;
    b = randn(n,1);
    
    %prepare the optimization: Constraints are x>0
    fctn = @(x) xSq(x,b);
    x0   = rand(n,1);
    xOpt = GaussNewtonProjected(fctn,x0,...
           'T1Min',0,'M0Min',0,'M0Max',Inf,'T1Max',Inf);
    
    %output:
    fprintf('\n\nSolved the problem\n\n\tMinimize: \n\t\tJ(x):=(x-b)''*(x-b)')
    fprintf('\n\t\t s.t. x>=0\nFor b=');
    b'
    fprintf('Solution was xOpt=');
    xOpt'
    
    
end



function [f,para,df,H] = xSq(x,b)

    f    = 1/2*((x-b)'*(x-b));
    df   = (x-b)';
    H    = eye(numel(x));
    para = struct();
end