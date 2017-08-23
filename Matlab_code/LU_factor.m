function [L,U]=LU_factor(A)
  % Here, factorization A=LU is done without pivoting,
  % permutations of the rows and columns in A.
  % This function overwrites L and U on A.
  
  % input A is an n by n matrix.
  
  % Output L is a unit lower triangular matrix.
  % Output U is non-singular upper matrix.
  
  %% Pre-defining matrices and indices
  s=size(A);
  n=s(1);
  
  %% Compute L and U, ovewrwrite on A.
  
  for i=1:(n-1) % need to do n-1 operations on A
    
    for j=(i+1):n
      A(j,i)=A(j,i)/A(i,i);
    end
    
    for j=(i+1):n
      for k=(i+1):n
        A(j,k)=A(j,k)-A(j,i)*A(i,k);
      end
    end
  end
  
  %% Construct L, copy values from A
  L=eye(n); % pre-define as identity matrix.
  for i=2:n
    for j=1:(i-1)
      L(i,j)=A(i,j);
    end
  end
  
  %% Construct U, copy values from A
  U=zeros(n); % pre-define as zero matrix.
  for i=1:n
    for j=i:n
      U(i,j)=A(i,j);
    end
  end
  
end
