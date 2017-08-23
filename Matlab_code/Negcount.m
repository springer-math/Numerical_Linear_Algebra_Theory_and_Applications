% ----------------------------------------
%Compute number of eigenvalues of a tridiagonal matrix A
%(without pivoting)  which are less then z
% ----------------------------------------

function [ neg ] = Negcount( A,z )
  
  d=zeros(length(A),1);
  d(1)=A(1,1)-z;
  for i = 2:length(A)
    d(i)=(A(i,i)-z)-(A(i,i-1)^2)/d(i-1);
  end
  
  %compute number of negative eigenvalues of A
  neg=0;
  for i = 1:length(A)
    if d(i)<0
      neg = neg+1;
    end
  end
  
end

