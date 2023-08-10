function [k,lamda_start] = GradSelectRule(Ndim)
% An evaluation rule for gradient selection
% -----------------------------------------
% input:
% Ndim: number of design variables

% output:
% k: number of active variables

if nargin==0
%     for test purpose
    Ndim=20;
end
e_exp1 = 1/(exp(1)-1);
part1 = 1/((1/(1+Ndim))^3-1);
F=@(lamda)log(lamda+e_exp1)-log(e_exp1)-(-part1+part1*((1+lamda*Ndim)/(1+Ndim))^3);
options = optimoptions('fsolve','Display','off');
lamda_start = fsolve(F,0.5,options);
k = round(lamda_start*Ndim);
end