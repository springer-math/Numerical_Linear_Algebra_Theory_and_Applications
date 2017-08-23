% ----------------------------------------
%   Solution of the system of linear equations Ax =  b via
%   SVD decomposition of a matrix A.
%   SVD decomposition is done via matlab function svd.
%   Matrix A is m-by-n, m > n, the vector of the rhs b is of the size n.
% ----------------------------------------

function x=LLSSVD(A,b)
  
  [U, S, V]=svd(A);
  
  UTb=U'*b;
  
  % choose tolerance
  tol=max(size(A))*eps(S(1,1));
  s=diag(S);
  n=length(A(1,:));
  
  % compute number of singular values > tol
  r=sum(s > tol);
  
  w=[(UTb(1:r)./s(1:r))' zeros(1,n-r)]';
  
  x=V*w;
