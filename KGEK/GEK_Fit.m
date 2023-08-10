function gekmodel = ...
    GEK_Fit(xi,y,erry,grad,errgrad,index_grad,options)
%% Fit the GEK model
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% xi - location of observations, size N x ndim
% y - observed QoI, size N x 1
% erry - measurement error in y, size N x 1
% grad - observed gradients, size N* x ndim*
% errgrad - measurement error in grad, size N* x ndim*
% index_grad - location of gradient exists, 0 or 1 matrix in size n x ndim
% options - GEK options, refer to defaultopts
%
% NOTE: specify '' for grad and errgrad if no gradient information is 
% available
%
% == OUTPUT ==
% gekmodel - struct containing model info

%% === add default options ===
options = defaultopts(options,xi);

if ~isempty(index_grad)
    index_grad = index_grad>0;
end
%% === normalize ===
[xin,yn,errn,norming,normreport] = ...
    normalize(xi,y,erry,grad,errgrad,index_grad,options);

%% === estimate hyperparameters ===
% this is still incomplete
[hypers,fit,options] = ...
    hyperestimate(xin,yn,errn,index_grad,options);

%% compile output
if isfield(options,'Ncomp')
    theta = hypers(3:2+options.Ncomp)';
else
    theta = hypers(3:end);
end

gekmodel = struct('regr',options.regression, 'theta',theta.', 'hypers',hypers, ...
  'beta',fit.beta, 'gamma',fit.gamma, 'sigma2',norming.Ysc(2,:).^2.*fit.sigma2, ...
  'S',xin,'Y',yn, 'Ssc',norming.Ssc, 'Ysc',norming.Ysc, 'F',fit.F,'R',fit.R, ...
  'C',fit.C, 'Ft',fit.Ft, 'G',fit.G,'index_grad',index_grad ,...
  'options',options ,'N',norming.N,'ModelType',normreport.ModelType);        

if isfield(fit,'F_del_index')
    gekmodel.F_del_index = fit.F_del_index;
end
if strcmp(options.hyperest,'none')% for WGEK
    gekmodel.beta1 = (fit.F)'*(fit.R\fit.F);
    gekmodel.beta2 = (fit.F)'*(fit.R\yn);
end
end
