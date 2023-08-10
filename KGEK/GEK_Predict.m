function [y, mse, dy, dmse] = GEK_Predict(gekmodel,x)
% Kriging/GEK predictor

%% prepare data
x           = (x - gekmodel.Ssc(1,:))./gekmodel.Ssc(2,:);
xi          = gekmodel.S;
theta       = gekmodel.theta; 
index_grad  = gekmodel.index_grad;
options     = gekmodel.options;
[m, ndim]   = size(xi);
if isfield(options,'Ncomp')
    % for KPLS
    coef_pls = options.coef_pls;
    theta = sum(theta.*coef_pls.^2,2)';
%     theta = (2*thetaKPLS).^(-0.5);%
end
%% calculate r and f
% calculate and delete unneccessary items in f for PGEK
[f,df] = feval(gekmodel.regr, x);
switch gekmodel.ModelType
    case 'Kriging'
        covars = 'no';
    otherwise %{'GEK','GGEK','PGEK','HGEK'}
        covars = 'yes';
end
% calculate r
rT = corrmatrix(x,xi,theta,options.predicttype,covars,options,index_grad);
rT = rT';
%% prediction and mse
% ------------------------------------
% Scaled predictor
sy = f * gekmodel.beta + (gekmodel.gamma * rT)';
%  Predictor
y = gekmodel.Ysc(1,:) + gekmodel.Ysc(2,:) * sy;
%----------------------------mean square error
rt = gekmodel.C \ rT;
u = gekmodel.G \ (gekmodel.Ft.' * rt - f.');
mse = gekmodel.sigma2 .* (1 + sum(u.^2,1) - sum(rt.^2,1))';
% mse = gekmodel.sigma2 .* (1 + u.^2 - sum(rt.^2,1))';

%% gradient of y and mse for GEK with zero-order regression only
if nargout>2 && strcmp(gekmodel.regr,'regpoly0')
        % implementation based on 
        % Optimization of expensive black-box problems via Gradient-enhanced Kriging
        mx = size(x,1);    % number of trial sites
        F = gekmodel.F;
        b = corrmatrix(x,xi,theta,'predictplus',covars,options,index_grad);
        b(1:mx,:) = [];% b(1:mx,:) is for prediction and mse, hence deleted
        % dy 
        % for regpoly0, df=0 
%         youtn = df * gekmodel.beta + b*(gekmodel.R\(gekmodel.Y));%correct
%         youtn = df * gekmodel.beta + b*(gekmodel.R\(gekmodel.Y-F*gekmodel.beta));
        youtn = df * gekmodel.beta + b*gekmodel.gamma';% based on DACE
        youtn = reshape(youtn,[],ndim);
        dy = gekmodel.Ysc(2,:)*youtn./gekmodel.Ssc(2,:);
        % dmse && strcmp(gekmodel.regr,'regpoly0')
        if nargout>3 
            maxrows = floor( options.maxArraySize / mx / mx );
            if mx > maxrows
                warning('max memory reached! loop to calculate dmse~')
                nloops = ceil(mx/maxrows);
                dmse = zeros(mx,ndim);
                for i = 1:nloops
                    from = (i-1)*maxrows+1;
                    to = min(i*maxrows,mx);
                    xtemp = x(from:to,:);
                    mxT = size(xtemp,1); 
                    b = corrmatrix(xtemp,xi,theta,'predictplus',covars,options,index_grad);
                    b(1:mxT,:) = [];
                    switch gekmodel.corr% corrgauss matern32 matern52
                        case 'corrgauss'
                            diagp = ones(mxT,1)*2*theta';
                        case 'matern52'
                            diagp = ones(mxT,1)*5/3*theta';
                    end
                    diagp = reshape(diagp,numel(diagp),1);
                    omit_term = ((F'*(gekmodel.R\b')).^2/(F'*(gekmodel.R\F)))';
                    varxoutn =  diagp - diag(b*(gekmodel.R\b')) + omit_term;
                    varxoutn = reshape(varxoutn,[],ndim);
                    dmseT = gekmodel.sigma2*varxoutn./gekmodel.Ssc(2,:).^2;
                    dmse(from:to,:) = dmseT;
                end
            else
                switch gekmodel.corr% corrgauss matern32 matern52
                    case 'corrgauss'
                        diagp = ones(mx,1)*2*theta';
                    case 'matern52'
                        diagp = ones(mx,1)*5/3*theta';
                end
                diagp = reshape(diagp,numel(diagp),1);
                omit_term = ((F'*(gekmodel.R\b')).^2/(F'*(gekmodel.R\F)))';
                varxoutn =  diagp - diag(b*(gekmodel.R\b')) + omit_term;
                varxoutn = reshape(varxoutn,[],ndim);
                dmse = gekmodel.sigma2*varxoutn./gekmodel.Ssc(2,:).^2;
            end
%         else
%             disp(' dmse calculation only support regpoly0 ')
        end
end
end
