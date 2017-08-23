% ----------------------------------------
%   Solution of the system of linear equations A^T Ax = A^T b
%   using Cholesky factorization of A^T A.
%   Matrix A is m-by-n, m > n, the vector of the rhs b is of the size n.
% ----------------------------------------

function x=LLSChol(A,b)
  
  ATb=A'*b;
  ATA=A'*A;
  n=length(A(1,:));
  lowerChol=zeros(n);
  
  %Cholesky factorization
  for j=1:1:n
    s1=0;
    for k=1:1:j-1
      s1=s1+lowerChol(j,k)*lowerChol(j,k);
    end
    lowerChol(j,j)=(ATA(j,j)-s1)^(1/2);
    for i=j+1:1:n
      s2=0;
      for k=1:1:j-1
        s2=s2+lowerChol(i,k)*lowerChol(j,k);
      end
      lowerChol(i,j)=(ATA(i,j)-s2)/lowerChol(j,j);
    end
  end
  
  % Solver for LL^T x = A^Tb:
  % Define z=L^Tx, then solve
  % Lz=A^T b to find z.
  % After by known z we get x.
  
  % forward substitution Lz=A^T b to obtain z
  
  for i=1:1:n
    for k=1:1:i-1
      ATb(i)=ATb(i)-ATb(k)*lowerChol(i,k);
    end
    ATb(i)=ATb(i)/lowerChol(i,i);
  end
  
  % Solution of L^Tx=z , backward substitution
  
  for i=n:-1:1
    for k=n:-1:i+1
      ATb(i)=ATb(i)-ATb(k)*lowerChol(k,i);
    end
    ATb(i)=ATb(i)/lowerChol(i,i);
  end
  
  % Obtained solution
  x=ATb;
  
