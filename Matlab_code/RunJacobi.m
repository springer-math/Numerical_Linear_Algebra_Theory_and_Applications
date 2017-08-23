% ----------------------------------------
% Run Classical Jacobi rotation algorithm.
% until the matrix A is sufficiently diagonal or off(A) < tol
% ----------------------------------------

function [A] = RunJacobi(A)
  
  tol=0.005;
  
  iter=1;
  
  %compute initial off's
  [sum,v]=off(A);
  
  while sum >tol && iter<100000
    % search for maximal values of off's
    j=v(2,max(v(1,:)) == v(1,:)); %get index j
    k=v(3,max(v(1,:)) == v(1,:)); %get index k
    
    %perform Jacobi rotation for indices (j,k)
    A=jacobiRot(A,j,k);
    [sum,v]=off(A);
    iter=iter+1;
  end
  
end

% Run one Jacobi rotation

function [A] = jacobiRot( A,j,k )
  tol=0.0001;
  
  if abs(A(j,k))>tol
    tau=(A(j,j)-A(k,k))/(2*A(j,k));
    t=sign(tau)/(abs(tau)+sqrt(1+tau^2));
    c=1/(sqrt(1+t^2));
    s=c*t;
    
    R=eye(length(A));
    R(j,j)=c;
    R(k,k)=c;
    R(j,k)=-s;
    R(k,j)=s;
    
    A=R'*A*R;
  end
  
end

% Compute  off's:  the square root of the sum of squares
% of the upper off-diagonal elements.
% v is a matrix that holds the information needed.

function [sum,v] = off(A)
  
  sum=0;
  %create array v for off's:
  % in the first row  will be sum of square root of the squares of computed off's
  % in the second row: the index  j
  % in the third row: the index k
  
  v=[0;0;0];
  for i=1:(length(A)-1)
    for j=(i+1):length(A)
      sum=sum+A(i,j)*A(i,j);
      v=[v,[sqrt(A(i,j)*A(i,j));i;j]];
    end
  end
  sum=sqrt(sum);
  v=v(:,2:end);
  
end

