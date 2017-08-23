% ----------------------------------------
%   Construction of columns in matrix A  using linear splines.
%   Input arguments: T - column vector with junction points,
%   x are measurement ponts (discretization points).
%   Returns column with number k to the matrix A.
% ----------------------------------------

function f=fihatt(k,x,T)
  
  h=diff(T);
  N=length(T);
  f=zeros(size(x));
  if k>1
    I=find(x>=T(k-1) & x<=T(k));
    f(I)=(x(I)-T(k-1))/h(k-1);
  end
  if k<N
    I=find(x>=T(k) & x<=T(k+1));
    f(I)=(T(k+1)-x(I))/h(k);
  end
  
