% ----------------------------------------
%   Solution of the system of linear equations Ax =  b via
%  QR decomposition of a matrix A.
%   Matrix A is m-by-n, m > n, the vector of the rhs b is of the size n.
%  QR decomposition of A is done via classical
%  Gram-Schmidt  (CGM) orthogonalization procedure.
% ----------------------------------------

function x=LLSQR(A,b)
  
  n=length(A(1,:));
  q=[];
  r=[];
  
  for i=1:1:n
    q(:,i)=A(:,i);
    for j=1:1:i-1
      r(j,i)=q(:,j)'*A(:,i);
      q(:,i)=q(:,i)-r(j,i)*q(:,j);
    end
    r(i,i)=norm(q(:,i));
    q(:,i)=q(:,i)/r(i,i);
  end
  
  % compute right hand side in the equation
  Rx=q'*b;
  
  % compute solution via backward substitution
  for i=n:-1:1
    for k=n:-1:i+1
      Rx(i)=Rx(i)-Rx(k)*r(i,k);
    end
    Rx(i)=Rx(i)/r(i,i);
  end
  
  x = Rx;
