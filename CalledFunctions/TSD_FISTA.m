%% Basic Information
% 
%	Input parameters:
%       A       : Measured system matrix
%       b       : Induced signal vector
%       d       : Subtracted noise vector (kappa*eta)
%       Phi     : Sum-difference frame/analysis operator (p*n)
%       mu      : Parameter of the total sum-difference regularization
% 
%   Output parameters:
%       x       : Desired concentration vector
% 
%   Parameters in opts: 
%       rtFlag	: Flag of real-time imaging (true) or not (false), the default is false
%       psFlag	: Flag of parameter searching (true) or not (false), the default is false
%       Frames	: Frames count, the default is 1
%       maxit	: Maximum number of iterations, the default is 1e3
%       maxt	: Maximum iteration time (s), the default is 0.1s (10 frames/second)
%       tol     : Tolerance of error, the default is 1e-5
%       zeta    : Algorithmic parameter, the default is 0.5
%       gamma   : Algorithmic parameter, the reciprocal of the smallest Lipschitz constant of nabla_F
%       T       : Temporary parameter for updating q quickly
%       f       : Function for reducing computational complexity
%       m       : Number of selected frequency components
%       n       : Number of pixels/voxels of reconstructed image
%       time1   : Time of data preprocessing
%       time2   : Time of iterations
%       iter    ：Final number of iterations
%       error   : Final relative error
%       record  : Record step parameters, x, err, satrt time of record, end time of record, index
% 
%   Email: jieh@buaa.edu.cn
% 
%   Copyright (C) May 2024 Jie He

%% Main Function: Projected General Forward-Backward Splitting (pGFBS) for the TSD Model
function [x,opts] = TSD_FISTA(A,b,d,Phi,mu,opts)
    
    %%%%% Operating Mode Selection and Algorithmic Parameter Definition %%%%%
    if ~isfield(opts, 'Frames'); opts.Frames = 1; end           % frames count
    if opts.Frames == 1
        if ~isfield(opts, 'rtFlag')                         	% flag of real-time imaging
            opts.rtFlag = false;                              	% the default is false
        elseif opts.rtFlag
            opts.psFlag = false;                               	% mutually exclusive
            opts.maxit = inf;                                 	% no iteration limits
        end       
        if ~isfield(opts, 'psFlag')                           	% flag of parameter searching
            opts.psFlag = false;                             	% the default is false
        elseif opts.psFlag
            opts.rtFlag = false;                               	% mutually exclusive
            opts.maxt = inf;                                  	% no time limits
        end  
        if ~isfield(opts, 'maxit'); opts.maxit = 1e3; end    	% maximum number of iterations
        if ~isfield(opts, 'maxt'); opts.maxt = 0.1; end       	% maximum iteration time (10 frames/second)
        if ~isfield(opts, 'tol'); opts.tol = 1e-5; end        	% tolerance of error
        if ~isfield(opts, 'SNR'); opts.SNR = 3; end          	% snr level
        if ~isfield(opts, 'zeta'); opts.zeta = 1e-3; end     	% an robust algorithmic parameter, we let it be 0.5            
    end
    
    %%%%% Data Preprocessing of the Main Function %%%%%
    tic;                                                        % start time of data preprocessing
    if opts.Frames == 1
        % 1) Frequency components selection %
        opts.freqs = abs(b) > (opts.SNR*abs(d));                % indexes of selected frequency components
        A = A(opts.freqs,:);                                    % final system matrix
        b = b(opts.freqs);                                      % final signal vector
        % 2) L2-norm calculation %
        Y1_2 = 1./vecnorm(A,2,2);                               % vectorization of row-weighting matrix of A
        A = A.*Y1_2;                                            % A=Y^{1/2}*A, row normalization of A
        [opts.m,opts.n] = size(A);                              % the size of system matrix
        if opts.m < opts.n  
            if isfield(opts, 'L')
            else
                opts.L = double(norm(single(A*A.'),2));             % norm(A*A.')=norm(A.'*A)
            end
        else
            G = A.'*A;                                          % nabla_F = G*x-A.'*Y^{1/2}*b, here A has been row-normalized
            if isfield(opts, 'L')
            else
                opts.L = double(norm(single(G),2));                    % smallest Lipschitz constant of nabla_F
                % opts.L = 2.1929e+03;
            end
        end
        % 3) Formula simplification %     
        opts.gamma = 1.0/(opts.L+1/opts.zeta);                  % 0<gamma<1/(L+1/opts.zeta)
        if opts.m < opts.n/3                                    % the ratio 3 is selected based on experience
            opts.f = @(v1,q) (A.'*(A*v1)+q);                    % computational complexity is O(mn)
        elseif opts.m > opts.n      
            opts.f = @(v1,q) (G*v1+q);                          % computational complexity is O(n^2)
        else
            G = A.'*A;                                          % parameter of function f (n*n)
            opts.f = @(v1,q) (G*v1+q);                         	% computational complexity is O(n^2)
        end 
        opts.T = A.'.*(opts.gamma.*transpose(Y1_2));           	% temporary parameter for updating q quickly
%         clearvars A Y1_2 G                                      % for saving memory, comment it out if you like
        q = opts.T*b;                                           % parameter of function f
        rho = opts.gamma*mu;                                    % for soft thresholding step
        Phi_T = Phi.';                                          % transpose of Phi
        opts.p = size(Phi,1);
        % 4) Augments initialization %
        x_hat = zeros(opts.p,1);                            	% initialization of x
        v_hat = zeros(opts.p,1);                             	% initialization of v
    else
        q = opts.T*b(opts.freqs);                               % update of q
        rho = opts.gamma*mu;                                    % update of rho
        Phi_T = Phi.';                                          % transpose of Phi
        x_hat = opts.x_hat;                                  	% initialization of x
        v_hat = opts.v_hat;                                  	% initialization of v
    end
    t1 = toc;                                                   % end time of data preprocessing

    %%%%% Iterations of the Main Function %%%%%  
    tic;                                                        % start time of iterations
    k = 0;                                                      % initialization of iterations k
    g = 1;
    v = Phi_T*v_hat;
    x = Phi_T*x_hat;
    while k<opts.maxit && toc<opts.maxt
        k = k+1;                                                % update of k  
        xp_hat = x_hat;
        xp = x;                                                 % record of x
        u_hat = (1-opts.gamma/opts.zeta)*v_hat+...
            Phi*(max((opts.gamma/opts.zeta)*v,0)+...
            opts.f(-opts.gamma*v,q));                           % update of x_hat, let eta=1                               
        x_hat = shrink(u_hat,rho);                              % update of x_hat, let eta=1
        gp = g;
        g = 0.5+sqrt(0.25+gp^2);
        gt = (gp-1)/g;
        v_hat = x_hat+gt*(x_hat-xp_hat);
        v = Phi_T*v_hat;
        x = Phi_T*x_hat;
        err = norm(xp-x)/(norm(xp)+1e-3);                       % update of error
        if err<opts.tol; break; end                             % stopping criterion
    end 
    x = max(x,0);                                               % projection onto nonnegative space 
    t2 = toc;                                                   % end time of iterations
    
    %%%%% Parameters Record %%%%% 
    opts.time1(opts.Frames) = t1;                               % record of t1
    opts.time2(opts.Frames) = t2;                               % record of t2
    opts.iterations(opts.Frames) = k;                           % record of final iterations
    opts.error(opts.Frames) = err;                              % record of final error
    if opts.rtFlag || opts.psFlag
        opts.Frames = opts.Frames+1;                            % update of frames
        if opts.rtFlag
            opts.x_hat = x_hat;                              	% record of x     
            opts.v_hat = v_hat;                                 % record of v
        else     
            opts.x_hat = zeros(opts.p,1);                     	% initialization of x
            opts.v_hat = zeros(opts.p,1);                    	% initialization of v
        end
    end    
end

%% Auxliary Functions
function y = shrink(x,threshold)
    y = sign(x).*max(abs(x)-threshold,0);
end

%% Projected General Forward-Backward Splitting (pGFBS) for TSD Model
% TSD model :
%
% $$\arg\min\limits_{x\in\mathbf{C}^{n}}\quad\frac{1}{2}\|Y^{\frac{1}{2}}(Ax-b)\|_{2}^{2} +
% \mathcal{I}_{+}(x) + \mu\|\Phi x\|_{1}.$$
%
% Here, $A\in\mathbf{C}^{m\times n}$ is the measured system matrix,
% $b\in\mathbf{C}^{m\times 1}$ the induced voltage vector,
% $x\in\mathbf{C}^{n\times 1}$ the desired concentration vector,
% $\mathcal{I}_{+}(x)$ the indicator function,
% $\Phi\in\mathbf{R}^{p\times n}$ the normalized weighted sum-difference
% frame, $\mu>0$ the weighting parameters of 
% total sum-difference regularization, respectively.
%
% Let
%
% $$\mathcal{F}(x)=\frac{1}{2}\|Y^{\frac{1}{2}}(Ax-b)\|_{2}^{2}.$$
%
% Then, the GFBS iteration reads as
%
% $$\nabla\mathcal{F}(x_{k})=A^{H}Y(Ax_{k}-b),$$
%
% $$u_{k+1}=u_{k}+\theta_{k}\left[\textnormal{prox}_{\frac{\gamma}{\zeta}\mathcal{I}_{+}}
% \Big(2x_{k}-u_{k}-\gamma\nabla\mathcal{F}(x_{k})\Big)-x_{k}\right]$$
%
% $$\cdots\thinspace\enspace=u_{k}+\theta_{k}\Big[\mathcal{P}_{+}\Big(2x_{k}-u_{k}-\gamma\nabla\mathcal{F}(x_{k})\Big)-x_{k}\Big],$$
%
% $$v_{k+1}=v_{k}+\theta_{k}\left[\textnormal{prox}_{\frac{\gamma\mu}{1-\zeta}\|\Phi\cdot\|_{1}}
% \Big(2x_{k}-v_{k}-\gamma\nabla\mathcal{F}(x_{k})\Big)-x_{k}\right]$$
% 
% $$\cdots\thinspace\enspace =v_{k}+\theta_{k}\left[\Phi^{H}\mathcal{T}_{\frac{\gamma\mu}{1-\zeta}}
% \Big(\Phi\big(2x_{k}-v_{k}-\gamma\nabla\mathcal{F}(x_{k})\big)\Big)-x_{k}\right],$$
%
% $$x_{k+1}=\zeta u_{k+1}+(1-\zeta)v_{k+1}.$$
%
% Here, $\gamma\in(0,\frac{2}{L})$, $\theta_{k}\in(0,\min(\frac{3}{2},\frac{1}{2}+\frac{1}{\gamma L}))$,
% $\zeta>0$ are the GFBS parameters. $\nabla$ denotes the gradient operator.
% $L$ is the smallest Lipschitz constant of $\nabla\mathcal{F}$,
% equals to $\|A^{H}YA\|_{2}$.
