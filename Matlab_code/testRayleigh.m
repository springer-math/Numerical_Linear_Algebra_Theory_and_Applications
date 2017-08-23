% ----------------------------------------
% Program which generates predefined random tridiagonal matrices A of dim(A)=n
% and then calls the function   RayleighQuotient.m
% ----------------------------------------

n=10;
A=zeros(n);

for i=2:n
  tal = rand*30;
  A(i,i)=rand*20;
  A(i,i-1)=tal;
  A(i-1,i)=tal;
end
A(1,1)=22*rand;

%run algorithm of Rayleigh Quotient Iteration

[rq]=RayleighQuotient(A);

disp('Computed  Rayleigh Quotient is:')
disp(rq)

