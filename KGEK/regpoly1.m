function  [f, df] = regpoly1(S)
%REGPOLY1  First order polynomial regression function
%
% Call:    f = regpoly1(S)
%          [f, df] = regpoly1(S)
%
% S : m*n matrix with design sites
% f = [1  s]
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update April 12, 2002

[m, n] = size(S);
f = [ones(m,1)  S];%2.22
if  nargout > 1
%   df = [zeros(n,1) eye(n)];
% added by HYW @2020-07-21
  df1 = [zeros(n,1) eye(n)];
  df=zeros(m*n,n+1);
  for i=1:n
      df((i-1)*m+1:i*m,:) = repmat(df1(i,:),m,1);
  end
end
end