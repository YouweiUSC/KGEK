function [state0,A,idnnan,initreport] = ...
    initialstate(xi,yn,errn,hypers,index_grad,options)
%% INITIALSTATE Initial state for GEK
% 
% [state0 A initreport] = ...
%    initialstate(xi,yn,errn,hypers,options)
%
% Computes the initial state: 
% 
%    state0 = A\yn, 
%
% where A = R + H'PH. Note that the initial state can contain value and
% gradient information, since yn is a (compiled) normalized data vector
%
% Calculation of the initial state is a expensive operation (deblurr or 
% pull-back). The initial state can later be used for cheap prediction
% of the QoI at xiout (blurr or push-forward or diffusion)
%
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% xi - location of observations, size N x ndim
% yn - normalized compiled data vector, size N x (1+ndim)
% errn - normalized error vector, size N x (1+ndim)
% hypers - optimized hyperparameters, size 1 x (ndim + 2)
% options - GEK options, refer to defaultopts
%
% == OUTPUT ==
% state0 - initial state, size N x (1+ndim)
% A - gain matrix, size N(1+ndim) x N(1+ndim)
% initreport - << under construction >>


%% == check for covariables ==
tic
[n1 m1] = size(xi);
[n2 m2] = size(yn);
if n1 == n2
    covars = 'no';
else
    covars = 'yes';
end

% if options.debug > 0
%     fprintf('gek -> finding initial state ...\n')
% end

%% == split hypers ==
errnf = nan(n2,1);
if n1 == n2
    errnf = hypers(1) * errn;
else
    errnf(1:n1) = hypers(1) * errn(1:n1);
    errnf(n1+1:end) = hypers(2) * errn(n1+1:end);
end

range = hypers(3:end);

%% == determine initial state ==

% build matrices
P = corrmatrix(xi,xi,range,'analyze',covars,options,index_grad);    % prior
R = diag(errnf.^2);                   % likelihood
A = R + P;                           % gain (is output for later use)

% remove nan's
idnnan = find(~isnan(yn));
yn = yn(idnnan);
A  = A(idnnan,idnnan);

% get initial state
state0 = A\yn;                       % cholesky -> backward substitution

%% == check the accuracy of the initial state ==
prec = eps;
accur = std(A*state0-yn);
if options.debug > 0
%     fprintf(['   residual : ' num2str(accur) '\n'])
    if accur > sqrt(prec);
        warning(['the initial state for size(xi) = [ ' num2str(size(xi)) ' ] is quite fuzzy'])
    end
end

%% == report ==
initreport.cpu = toc;
initreport.accuracy = accur;

if options.debug > 2
    if strcmp(covars,'no'), figure(11), else figure(11), end
    subplot(2,1,1),contourf(A,100,'edgecolor','none'), colorbar
    if strcmp(covars,'no'), title('kriging A'), else title('gek A'), end
    subplot(2,1,2),contourf(chol(A),100,'edgecolor','none'), colorbar
    if strcmp(covars,'no'), title('kriging chol(A)'), else title('gek chol(A)'), end
end
