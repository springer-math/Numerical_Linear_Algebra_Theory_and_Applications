% ----------------------------------------
% generation of the random  tridiagonal symmetric matrix
% ----------------------------------------

function [A] = randomTridiag(n)
  
  A=zeros(n);
  
  for i=2:n
    num = rand*30;
    A(i,i)=rand*20;
    A(i,i-1)=num;
    A(i-1,i)=num;
  end
  A(1,1)=22*rand;
end

% ----------------------------------------
