% ----------------------------------------
% Computes value of  Rayleigh Quotient rq which is in the tolerance
% tol from an eigenvalue of A
% ----------------------------------------

function rq = RayleighQuotient(A)
  
  [n,~]=size(A);
  x0=zeros(n,1);
  
  % initialize initial vector x0 which  has norm 1
  x0(n)=1;
  
  tol = 1e-10;
  xi = x0/norm(x0,2);
  
  i=0;
  % initialize  Rayleigh Quotient for x0
  rq = (xi'*A*xi)/(xi'*xi);
  
  while norm((A*xi-rq*xi),2) > tol
    yi = (A-rq*eye(size(A)))\xi;
    xi=yi/norm(yi,2);
    rq = (xi'*A*xi)/(xi'*xi)
    i=i+1;
  end
  
end

% ----------------------------------------
