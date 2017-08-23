% ----------------------------------------
%                Inverse Iteration or Inverse Power method
% Computes eigenvalue closest to sigma  and corresponding eigenvector
% ----------------------------------------
clc
clear all
close all
eps = 1e-17;
fig = figure;

for i =1:4
  if(i==1)
    % Matrix not diagonalizable
    n=2;
    A =[0 10;0 0];
    % Matrix has two real eigenvalues with the same sign
    % n=3;
    % A =[5 0 0;0 2 0;0 0 -5];
  elseif (i==2)
    % Matrix has four real eigenvalues with the same sign
    n =4;
    A=[3,7,8,9;5,-7,4,-7;1,-1,1,-1;9,3,2,5];
  elseif (i ==3)
    % Largest eigenvalue is complex
    n =3;
    A =[0 -5 2; 6 0 -12; 1 3 0];
  elseif (i==4)
    % n =2;
    % A =[7 -2;3 0];
    n=5;
    A=rand(5,5);
  end
  
  % get  reference values of eigenvalues
  exact_lambda = eig(A);
  
  %make orthogonalization
  Q=orth(rand(n,n));
  
  A= Q'*A*Q;
  
  % set initial guess for the eigenvector x0
  x0=rand(n,1);
  x0=x0/norm(x0);
  lambda0 = inf;
  % choose a shift:  should be choosen as closest to the desired eigenvalue
  sigma=10;
  % lambda1 = 0;
  lambdavec =[];
  count =1;
  % main loop in the power method
  while (count <1000)
    A_shift = A - sigma*eye(size(A));
    y1= inv(A_shift)*x0;
    x1=y1/norm(y1);
    lambda1 = transpose(x1)*A*x1;
    lambdavec(count)= lambda1 ;
    x0=x1;
    if(abs(lambda1 - lambda0 )<eps )
      break ;
    end
    lambda0 = lambda1 ;
    count = count + 1;
  end
  
  % Print  computed and exact  eigenvalue
  str =['Example ' num2str(i)];
  str =[str, '. Comp. eig.:' num2str(lambda1)];
  str=[str, ', Ref. eig.:' num2str(exact_lambda',2)];
  subplot (2 ,2,i)
  plot (lambdavec,  'LineWidth',3)
  xlabel(' Number of iterations in Inverse Power method')
  ylabel('Computed eigenvalues')
  title (str, 'fontsize',12)
  
end
