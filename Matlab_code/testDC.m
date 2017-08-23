% ----------------------------------------
% Program which generates predefined random tridiagonal matrices A of dim(A)=n
% and then calls the function   DivideandConq.m
% ----------------------------------------

%Program which generates some random  symmetric tridiagonal matrices

n=5;
A=zeros(n);

for i=2:n
  tal = rand*30;
  A(i,i)=rand*20;
  A(i,i-1)=tal;
  A(i-1,i)=tal;
end
A(1,1)=22*rand;

%run Divide-and-Conquer  algorithm
[Q,L]=DivideandConq(A)

