function [xout, varxout, gradout, vargradout, denormreport] = ...
    denormalize(xiout,xoutn,varxoutn,norming,options)
%% DENORMALIZE De-normalize data for GEK
% 
% [xout varxout gradout vargradout denormreport] = ...
%    denormalize(xiout,xoutn,varxoutn,norming,options)
%
% De-normalizes data for GEK. In case of gradients, de-compiles value and
% gradient information individual vectors.
%
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% xiout - location of predictions, size n x ndim
% xoutn - predicted normalized mean, size n x 1
% varxoutn - predicted normalized variance, size n x 1
% norming - struct containing normalization info
% options - GEK options, refer to defaultopts
%
% == OUTPUT ==
% xout - predicted mean of QoI, size n x 1
% varxout - predicted variance of QoI, size n x 1
% gradout - predicted gradients, size n x ndim
% vargradout - predicted variance of gradients, size n x ndim
% denormreport - << under construction >>

tic;
%% == general ==
[nx, ndim] = size(xiout);

% if options.debug > 0
%     fprintf('gek -> denormalizing ...\n')
% end

%% == denormalize ==
% mean
ylin = norming.c(1) + xiout*norming.c(2:end);
xout = norming.st*xoutn(1:nx) + ylin;

% variance
if strcmp(options.estvar,'yes')
    varxout = norming.st^2 * varxoutn(1:nx);
elseif strcmp(options.estvar,'cheap')
    varxout = norming.st^2 * varxoutn(1);
else
    varxout = '';
end

% gradients
if strcmp(options.predicttype,'predictplus')
    gradout = norming.st*xoutn(nx+1:end);
    gradout = reshape(gradout,nx,ndim);
    vargradout = norming.st^2*varxoutn(nx+1:end);
    vargradout = reshape(vargradout,nx,ndim);
else
    gradout = NaN;
    vargradout = NaN;
end

%% == report ==
denormreport.cpu = toc;
