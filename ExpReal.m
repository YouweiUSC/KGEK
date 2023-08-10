clear, close all;clc;
addpath('KGEK','Ddcor','Results');% ,'WGEK','DACE'
% test and validation of GEK and GGEK's prediction and mse's derivative
%%
n_run = 1;   
VariableNamesT = {'Kriging','GEK','KGEK'};
N_method = length(VariableNamesT);
R2 = zeros(n_run,N_method);RMSE = R2; MAE = R2; ElaspeTime = R2;
%% settings
for icase=1
Prob = 'NACA0012';
ndim = 18;
sam_data = 'NACA0012_18.mat';
%%
disp(' overall step 0: generate test samples ')
%% --------------------------------------
for irun=1:n_run
    % generate samples 
    load(sam_data, 'sam_x', 'sam_y', 'sam_grad')
    N_split = 40;
    index = randperm(size(sam_y,1),N_split);
    test_x = sam_x;
    test_y = sam_y;
    sam_x = sam_x(index,:);
    sam_y = sam_y(index,:);
    sam_grad = sam_grad(index,:);
    test_x(index,:) = [];
    test_y(index,:) = [];
    disp(' step 1: generate samples ')
    N = size(sam_x,1); 
    xi = sam_x;
    y = sam_y(:,1);
    %% --------------------------------------
    % calculate gradients
    disp(' step 2: calculate gradients ')
    grad = sam_grad;
    % options for Kriging and GEK
    erry = 0;errgr = 0;
    options.hyperest = 'ga'; 
    %% --------------------------------------
    % model construction
        disp(' step 3: model construction ')
        model = cell(N_method,1);
        disp('kriging')
        tic%kriging
        model{1} = GEK_Fit(xi,y,erry*ones(size(y)),'','','',options);
        ElaspeTime(irun,1) = toc;
        disp('GEK')
        tic%GEK
        model{2} = GEK_Fit(xi,y,erry*ones(size(y)),grad,errgr*ones(size(grad)), ...
            ones(N,ndim),options);
        ElaspeTime(irun,2) = toc;
        disp('KGEK')
        %KGEK
        model{3} = KGEK_Fit(xi,y,grad,options);
        ElaspeTime(irun,3) =  model{3}.HyperEstTime;
    %% --------------------------------------
    % prediction
    disp(' step 4: prediction ')
    y_pred = zeros(size(test_y,1),N_method);
    y_pred(:,1) = GEK_Predict(model{1},test_x);
    y_pred(:,2) = GEK_Predict(model{2},test_x);
    y_pred(:,3) = KGEK_Predict(model{3},test_x);
    %% --------------------------------------
    % accuracy metric calculation
    disp(' step 5: accuracy metric calculation ')
    for iM = 1:N_method
        [R2(irun,iM),RMSE(irun,iM),MAE(irun,iM)] = ...
            ModelAccuracyMetric(y_pred(:,iM),test_y(:,1));
    end
end
%% metric statistic
RMSE(n_run+1,:) = mean(RMSE(1:n_run,:),1);
MAE(n_run+1,:)  = mean(MAE(1:n_run,:),1);
R2(n_run+1,:)   = mean(R2(1:n_run,:),1);
ElaspeTime(n_run+1,:)  = mean(ElaspeTime(1:n_run,:),1);
RMSE(n_run+2,:) = std(RMSE(1:n_run,:),1);
MAE(n_run+2,:)  = std(MAE(1:n_run,:),1);
R2(n_run+2,:)   = std(R2(1:n_run,:),1);
ElaspeTime(n_run+2,:)  = std(ElaspeTime(1:n_run,:),1);
rownamesT=cell(n_run+2,1);
for iname = 1:n_run
    rownamesT{iname} = ['run' num2str(iname)];
end
rownamesT{n_run+1} = 'mean';
rownamesT{n_run+2} = 'std';
R2t = array2table(R2,'rownames',rownamesT, 'VariableNames',VariableNamesT);
RMSEt = array2table(RMSE,'rownames',rownamesT, 'VariableNames',VariableNamesT);
MAEt = array2table(MAE,'rownames',rownamesT, 'VariableNames',VariableNamesT);
Timet = array2table(ElaspeTime,'rownames',rownamesT, 'VariableNames',VariableNamesT);
end