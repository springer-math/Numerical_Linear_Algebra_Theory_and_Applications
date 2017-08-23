% ----------------------------------------
% Computes the SVD decomposition of the matrix G
%  using the one-sided Jacobi rotation.
% ----------------------------------------

function [U,S,V] = RunSVDJacobi(G)
  
  % input tolerance
  tol=0.005;
  
  J=eye(length(G));
  iter=1;
  
  [sum,v]=off(G'*G);
  
  while sum>tol && iter<1000
    for j=1:(length(G)-1)
      for k=j+1:length(G)
        [G,J]=oneSidedJacobiRot(G,J,j,k);
      end
    end
    
    [sum,v]=off(G'*G);
    iter=iter+1;
  end
  
  % elements in the matrix sigma will be the two-norm
  % of i-column of the matrix G
  
  for i=1:length(G)
    sigma(i)=norm(G(:,i));
  end
  
  U=[];
  
  for i=1:length(G)
    U=[U,G(:,i)/sigma(i)];
  end
  
  V=J;
  
  S=diag(sigma);
  
end

% compute one-sided Jacobi rotation for G

function [G,J] = oneSidedJacobiRot(G,J,j,k )
  
  tol=0.0001;
  A=(G'*G);
  ajj=A(j,j);
  ajk=A(j,k);
  akk=A(k,k);
  
  if abs(ajk)>tol
    tau=(ajj-akk)/(2*ajk);
    t=sign(tau)/(abs(tau)+sqrt(1+tau^2));
    c=1/(sqrt(1+t^2));
    s=c*t;
    
    R=eye(length(G));
    R(j,j)=c;
    R(k,k)=c;
    R(j,k)=-s;
    R(k,j)=s;
    
    G=G*R;
    
    % if eigenvectors are desired
    J=J*R;
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
