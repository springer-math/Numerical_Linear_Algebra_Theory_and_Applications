function L=Cholesky(A)
  % Function factorizes square matrix A, assuming that A is s.p.d. matrix,
  % into A=LL', where L' is the transpose
  % of L, and L is non-singular lower triangular matrix.
  
  %%
  s=size(A);
  n=s(1);
  L=zeros(n);
  
  % diagonal elements i=j
  % a_jj=v_j*v_j'=l_j1^2+l_j2^2+...+l_jj^2 (sum has j-terms)
  
  % elements below diagonal, i>j
  % a_ij=v_i*v_j'=l_i1 l_j1 + l_i2 l_j2 + ... + l_ij l_jj (sum has j terms)
  
  for j=1:n % go through column 1 to n
    % Compute diagonal elements, i=j
    L(j,j)=A(j,j);
    for k=1:(j-1)
      L(j,j)=L(j,j)-L(j,k)^2;
    end
    L(j,j)=L(j,j)^(1/2);
    % Compute elements below diagonal, i>j
    for i=(j+1):n
      L(i,j)=A(i,j);
      for k=1:(j-1)
        L(i,j)=L(i,j)-L(i,k)*L(j,k);
      end
      L(i,j)=L(i,j)/L(j,j);
    end
  end
end
