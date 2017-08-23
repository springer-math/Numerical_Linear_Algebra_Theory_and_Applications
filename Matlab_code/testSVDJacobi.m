% ----------------------------------------
% Program which generates predefined random tridiagonal matrices A
% and the calls the function RunSVDJacobi.m
% ----------------------------------------

n=5;
A=zeros(n);

for i=2:n
  tal = rand*30;
  A(i,i)=rand*20;
  A(i,i-1)=tal;
  A(i-1,i)=tal;
end
A(1,1)=22*rand;

Ainit=A

disp('computed by one-sided Jacobi algorithm  SVD decomposition:');
[U,S,V]= RunSVDJacobi(Ainit)

disp('computed  SVD decomposition using svd command (for comparison):');
[u,sigma,v]=svd(Ainit)

