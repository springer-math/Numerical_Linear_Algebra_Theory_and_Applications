function x=ForwSub(L,b)
  % This function computes the vector $x$, of length $n$,
  % given $Lx=b$ where $L$ is an $n \times n$ nonsingular lower triangular matrix
  % and $b$ is a known vector of length $n$, by using forward substitution.
  
  %% Compute $x$ by forward substitution.
  s=size(L);
  n=s(1);
  x=zeros(n,1);
  % $L(i,i)*x(i)=b(i) - \sum_{j=1}^{i-1}$
  % First, set $x(i)=b(i)$, then subtract the known values.
  % Lastly, divide by diagonal entry $L(i,i)$
  x(1)=b(1)/L(1,1);
  for i=2:n
    x(i)=(b(i)-L(i,1:(i-1))*x(1:(i-1)))/L(i,i);
  end
end
