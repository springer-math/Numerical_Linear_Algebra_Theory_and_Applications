% ----------------------------------------
%                Power method
% ----------------------------------------
clc
clear all
close all
eps = 1e-7;
fig = figure;

for i =1:4
  if(i==1)
    % Matrix not diagonalizable
    % n=2;
    % A =[0 10;0 0];
    % Matrix has two real eigenvalues with the same sign
    n=3;
    A =[5 0 0;0 2 0;0 0 -5];
  elseif (i==2)
    % Matrix has four real eigenvalues with the same sign
    n =4;
    A=[3,7,8,9;5,-7,4,-7;1,-1,1,-1;9,3,2,5];
  elseif (i ==3)
    % Largest eigenvalue is complex
    n =3;
    A =[0 -5 2; 6 0 -12; 1 3 0];
  elseif (i==4)
    n =2;
    A =[7 -2;3 0];
    n=5;
    A=rand(n);
    
  end
  
  % get  reference values of eigenvalues
  exact_lambda = eig(A);
  
  % set initial guess for the eigenvector x0
  x0=rand(n,1);
  x0=x0/norm(x0);
  lambda0 = inf ;
  % lambda1 = 0;
  lambdavec =[];
  % counter for number of iterations
  count =1;
  % main loop in the power method
  
  while (count <1000)
    
    y1=A*x0;
    
    %  compute approximate eigenvector
    
    x1=y1/norm(y1);
    
    %  compute approximate eigenvalue
    lambda1 = transpose(x1)*A*x1;
    
    lambdavec(count)= lambda1 ;
    x0=x1;
    if(abs(lambda1 - lambda0 )<eps )
      break ;
    end
    lambda0 = lambda1;
    count = count + 1;
    
  end
  
  % Print  computed eigenvalue
  str =['Computed eigenvalue:' num2str(lambda1)];
  str=[str, ', Exact eigenvalues:' num2str(exact_lambda',2)];
  subplot (2 ,2,i)
  plot (lambdavec,  'LineWidth',3)
  xlabel('Number of iterations in Power method')
  ylabel('Computed eigenvalue')
  title (str, 'fontsize',10)
  
end
