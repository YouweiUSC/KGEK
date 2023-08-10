function options = ...
    defaultopts(options,xi)
%% DEFAULTOPTS Complements options for GEK
% 
% options = defaultopts(options,xi)
%
% Complements options with default GEK options. Please refer to code below
% for more details.
%
% == INPUT ==
% options - struct with options
% xi - location of observations, size N x ndim
%
% == OUTPUT ==
% options - struct with options
%

[n, ndim] = size(xi);

%% general
if isfield(options,'debug') == 0
    options.debug = 0; % 0 for no debug
                       % 1 for command line feedback
                       % 2 for cl feedback & plots
                       % -1 for test-phase code
end
% check specification
if isfield(options,'Ncomp')
    if options.Ncomp>ndim
    error('No. of component for KPLS must not be larger than No. of Variable')
    end
end
%% regression and correlation function
if isfield(options,'regression') == 0
    options.regression = 'regpoly0';% regpoly0 regpoly1 regpoly2
end

if isfield(options,'corr') == 0
    options.corr = 'matern32';% corrgauss matern32 matern52
end

%% hyperestimate.m
if isfield(options,'hyperinit') == 0
%     if ~isfield(options,'Ncomp')
        options.hyperinit = [1 1 ones(1,ndim)];
%     else
%         options.hyperinit = [1 1 ones(1,options.Ncomp)];
%     end
end

if isfield(options,'hyperspace') == 0
%     if ~isfield(options,'Ncomp')
        options.hyperspace = [0 0 ones(1,ndim)];
%     else
%         options.hyperspace = [0 0 ones(1,options.Ncomp)];
%     end
end

if isfield(options,'LB_hyper') == 0
    options.LB_hyper = 1e-6;
end

if isfield(options,'UB_hyper') == 0
    options.UB_hyper = 1e1;
end

if isfield(options,'hyperest') == 0
%     options.hyperest = 'none';
%     options.hyperest = 'fmin';
    options.hyperest = 'ga';
end

if isfield(options,'goalfun') == 0
%     options.goalfun = 'mle';
    options.goalfun = 'mle_robust';
end

if isfield(options,'fminopts') == 0
    options.fminopts = ...
        optimset('Display','on','TolFun',1e-20,...
            'TolX',1e-3,'MaxIter',1e5,'MaxFunEvals',2000,'Algorithm','sqp');
end

if isfield(options,'GAopts') == 0
    options.GAopts = optimoptions('ga','MaxGenerations',125,'MaxStallGenerations',50, ...
            'PopulationSize',8*ndim,'CrossoverFraction',0.8,'MigrationFraction',0.2, ...
            'Display','off');
end

if isfield(options,'optionsFminbnd') == 0
    options.optionsFminbnd = optimset('Display','on','MaxFunEvals',200,'TolX',1e-5);
end
%% predict.m
if isfield(options,'estvar') == 0
    options.estvar = 'yes';
end

if isfield(options,'predicttype') == 0
    options.predicttype = 'predict';% predictplus
end

if isfield(options,'maxArraySize') == 0
    if exist('memory') == 0
        options.maxArraySize = 1e8;
    else
       user = memory;
       safety = 1.3;
       options.maxArraySize = user.MaxPossibleArrayBytes / 8 / safety;
    end
end
