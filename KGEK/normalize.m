function [xi,yn,errn,norming,normreport] = ...
    normalize(xi,y,erry,grad,errgrad,index_grad,options)
%% NORMALIZE Normalize data for GEK
% 
% [yn errn norming N normreport] = ...
%    normalize(xi,y,erry,grad,errgrad,options)
%
% Normalizes data for GEK. In case of gradients, compiles value and
% gradient information to single data vector.
%
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
% == OUTPUT ==
% yn - normalized data vector, size N x (1+ndim)
% errn - normalized error vector, size N x (1+ndim)
% norming - struct containing normalization info
% normreport - << under construction >>

tic
%% == general ==
[ny, ndim] = size(xi);
N = ny;
if strcmp(grad,'')
    ModelType='Kriging';
elseif sum(sum(index_grad,1))<ny*ndim
    ModelType='GGEK';
elseif sum(sum(index_grad,1))==ny*ndim
    ModelType='GEK';    
end
if options.debug > 2
    fprintf('gek -> normalizing ...\n')
end

%% == normalize data ==
mS = mean(xi);   sS = std(xi);
% mS = 0;   sS = 1;
mY = mean(y);    sY = std(y);
% normalization
xi = (xi-mS)./sS;
ynn  = (y-mY)./sY;
% ynn = (y-ylin)/norming.st;
errynn = erry/sY;
if ~strcmp(grad,'')
    gradnn = sS.*grad/sY;
    errgradnn = sS.*errgrad/sY;
%     N = (ndim+1)*N;
    N = sum(sum(index_grad,1));
end
norming.N = N;
norming.Ssc = [mS; sS];
norming.Ysc = [mY; sY];
% compile data
if strcmp(grad,'')
    yn = ynn;
    errn = errynn;
else
    gradnn = reshape(gradnn(index_grad),N,1);
    errgradnn = reshape(errgradnn(index_grad),N,1);
    yn = [ynn ; gradnn];
    errn = [errynn ; errgradnn];
end

%% == report ==
normreport.cpu = toc;
normreport.ModelType = ModelType;
end
