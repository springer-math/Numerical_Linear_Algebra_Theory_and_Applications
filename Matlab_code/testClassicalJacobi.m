% ----------------------------------------
% Program which generates predefined random tridiagonal matrices A
% and the calls the function RunJacobi.m
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

% initialization of matrix
%A=rand(5,5)*10;

Ainit=A
%Ainit =A*A'

% run classical Jacobi algorithm
A= RunJacobi(Ainit)

%Print out computed by Jacobi algorithm eigenvalues
disp('computed by Jacobi algorithm eigenvalues:');
eig(A)

% Print out eigenvalues of the initial matrix A using eig(Ainit)
disp('eigenvalues of the initial matrix Ainit using eig(Ainit):');
eig(Ainit)

