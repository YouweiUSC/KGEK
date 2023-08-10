function PGEK = KGEK_Fit(xi,y,grad,options)
%% options
erry = 0;
errgr = 0;
[N,Ndim] = size(xi);
tic
%% calculate MIC
MIC = zeros(1,Ndim);
for idim = 1:Ndim
    MIC_sta = mine(xi(:,idim)', y');
    MIC(idim) = MIC_sta.mic;
end
%% specify options for auxiliary Kriging
optionsPreKrig = options;
optionsPreKrig.hyperest = 'mic';
optionsPreKrig.IxyW = MIC;
% optionsPreKrig.optionsFminbnd = optimset('Display','off','MaxFunEvals',1000,'TolX',1e-5);
disp('   build fast kriging')
Kriging = GEK_Fit(xi,y,erry*ones(size(y)),'','','',optionsPreKrig);
%% build GEK with full gradients
options.hyperinit = [1 1 reshape(Kriging.theta,1,Ndim)];
options.hyperest = 'fmin'; % ga fmin none
options.fminopts = ...
        optimset('Display','on','TolFun',1e-20,...
            'TolX',1e-3,'MaxIter',1e5,'MaxFunEvals',2000,'Algorithm','sqp');
disp('   build fast GEK')
PGEK = GEK_Fit(xi,y,erry*ones(N,1),grad,errgr*ones(size(grad)),ones(N,Ndim),options);
PGEK.MutualInformation = MIC;
PGEK.ModelType = 'KGEK';
% PGEK.var_active = I(1:k)';
% PGEK.N_active = k;
PGEK.HyperEstTime = toc;
end