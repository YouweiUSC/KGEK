function [R2,RMSE,MAE] = ModelAccuracyMetric(Y_hat,Y)
%
% Two input arguments: 'Y_hat' and 'Y'
% One output argument: 'R2'
%
% X:    Vector of x-parameter
% Y:    Vector of y-paramter
% R2:   Coefficient of determination

% Limitations
if length(Y_hat) ~= length(Y), error('Vector should be of same length');end
if nargin < 2, error('Not enough input parameters');end
if nargin > 2, error('Too many input parameters');end

N = length(Y_hat);
% R2: Coefficient of determination
ads_dif = abs(Y - Y_hat);
SSreg = sum(ads_dif.^2);
SStot = sum((Y - mean(Y)).^2);
R2 = 1-SSreg/SStot;

RMSE = sqrt(sum((Y-Y_hat).^2)/N);
MAE = max(abs(Y-Y_hat));
end