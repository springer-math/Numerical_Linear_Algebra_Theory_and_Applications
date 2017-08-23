function A=DiscretePoisson2D(n)
  % The function for 2D discretization  of the Laplace operator with sign minus: - laplace
  % Input parameters:
  % n -  number of inner nodes, which is assumed to be the same in both
  % the x_1- and x_2 directions.
  
  A = zeros(n*n,n*n);
  
  % Main diagonal
  for i=1:n*n
    A(i,i)=4;
  end
  
  % 1st and 2nd off-diagonals
  for k=1:n % go through block 1 to n
    for i=1:(n-1)
      A(n*(k-1)+i,n*(k-1)+i+1)=-1; %
      A(n*(k-1)+i+1,n*(k-1)+i)=-1;
    end
  end
  
  % 3rd and 4th off-diagonals
  for i=1:n*(n-1)
    A(i,i+n)=-1;
    A(i+n,i)=-1;
  end
  
end
