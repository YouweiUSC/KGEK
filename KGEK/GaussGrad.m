function  dr = GaussGrad(theta, d)
%CORRGAUSS  Gaussian correlation function,
%
%           n
%   r_i = prod exp(-d_ij^2/(2*theta_j ^2)) ,  i = 1,...,m
%          j=1
%
%
% Call:    dr = corrgauss(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 

[n, k] = size(d);  % number of differences and dimension of data

td = -0.5* d.^2 .* repmat(theta.^-2,n,1);
r = exp(sum(td, 2));

n_dr = n*(k+1);
dr = zeros(n_dr,k);
dr(1:n,:) = repmat(theta.^-2,n,1) .* d .* repmat(r,1,k);
for l=1:k
    for m = 1:k
        for j=1:n
            if m==l
                dr(m*n+j,l) = (theta(l)^-2-theta(l)^-4*d(j,l)^2)*r(j);
            else
                dr(m*n+j,l) = -theta(l)^-2.*theta(m)^-2*d(j,l)*d(j,m)*r(j);
            end
        end
    end
end