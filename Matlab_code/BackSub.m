function x=BackSub(U,b)
  % This function computes the vector $x$ by backward substitution.
  % We solve $Ux=b$, where $U$ is an $n \times n$ nonsingular upper triangular matrix
  % and $b$ is a known vector of the length $n$, finding the vector $x$.
  
  %% Compute x by backward substitution.
  s=size(U);
  n=s(1);
  x=zeros(n,1);
  %  $U(i,i)*x(i) = b(i) - \sum_{j=i+1}^{n}$
  x(n)=b(n)/U(n,n);
  for i=n-1:-1:1
    x(i)=(b(i)-U(i,(i+1):n)*x((i+1):n))/U(i,i);
  end
end
