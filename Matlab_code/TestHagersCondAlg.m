% ----------------------------------------
% Hager's algorithm: for the input matrix A
% the function HagerCond(A) computes
% the lower bound of the one-norm of the matrix A.
% ----------------------------------------

% First we generate some random  symmetric  matrices

n=5;
A=zeros(n);

for i=1:n
  for j=1:n
    tal = rand*30;
    A(i,i)=rand*20;
    A(i,j)=tal;
    A(j,i)=tal;
  end
end
disp(' The input matrix A  is:');

A

disp(' The computed lower bound of ||A||_1  is:');
HagersEst =  HagersAlg(A)

disp('  result of norm(A,1) is:');
norm(A,1)
