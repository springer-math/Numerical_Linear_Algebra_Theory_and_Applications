% ----------------------------------------
%  Matlab program HessenbergQR.m: we will first reduce the matrix A to the upper
%  Hessenberg matrix and then compute it's QR factorization
% ----------------------------------------
clc
clear all
close all
eps = 1e-07;
fig = figure;
N=10;
for i =1:6
  if(i==1)
    n=N;
    A=hilb(N);
  elseif (i==2)
    n=20;
    A=hilb(20);
  elseif (i ==3)
    % Largest eigenvalue is complex
    n =3;
    A =[0 -5 2; 6 0 -12; 1 3 0];
  elseif (i==4)
    % Matrix has four real eigenvalues
    n =4;
    A=[3,7,8,9;5,-7,4,-7;1,-1,1,-1;9,3,2,5];
  elseif (i==5)
    n =5;
    %
    A=[3,7,8,9,12;5,-7,4,-7,8;1,1,-1,1,-1;4,3,2,1,7;9,3,2,5,4];
  elseif (i==6)
    n=N;
    A= rand(N,N);
    
  end
  
  lambda0= inf(n,1);
  count = 1;
  iter =1;
  
  % get  exact eigenvalues in sorted order
  exact_lambda = sort(eig(A));
  
  % First we reduce matrix A to upper Hessenberg
  
  for k=1:n - 2
    x= A(k+1:n,k);
    u=x;
    u(1) = u(1)+ sign(x(1))*norm(x);
    u=u/norm (u);
    P= eye(n-k) - 2*(u*u') ;
    A(k +1:n,k:n) =P*A(k +1:n,k:n) ;
    A(1:n,k +1:n)=A(1:n,k+1:n)*P;
  end
  
% ----------------------------------------
  
  for k = 1:1000
    [Q,R] = qr(A);
    A = R*Q;
  end
  
  %%%%%%%%%% Find eigenvalues from Real Schur block
  j =2; count =1;
  eigs = zeros(1,n);
  while (j <=n)
    %real eigenvalues
    if(abs(A(j,j-1)) < 1e-3)
      eigs(j-1) =A(j -1,j -1);
      count= j -1;
    else
      % Complex  eigenvalues
      eigs(j-1: j)= eig(A(j -1:j,j -1:j));
      count =j;
      j=j +1;
    end
    j=j +1;
  end
  if(count < length(eigs))
    eigs(n)=A(n,n);
  end
% ----------------------------------------
  
  computed_lambda = sort(eigs)';
  
  if(norm(abs(computed_lambda - lambda0 ))<eps )
    break ;
  end
  lambda0 = computed_lambda ;
  iter = iter + 1;
  
  str =['Comp. eig.:' num2str(computed_lambda')];
  str=[str, ', Ex. eig.:' num2str(exact_lambda',2)];
  str_xlabel = ['Example ',num2str(i), '. Nr. of it. in method of QR it.:', num2str(iter)];
  
  subplot (3,2,i)
  
  plot (exact_lambda,'o b','LineWidth',2,'Markersize',10)
  hold on
  plot (computed_lambda,'+ r','LineWidth',2, 'Markersize',10)
  
  % xlabel(str, 'fontsize',10)
  
  xlabel('Real part of eigenvalues');
  ylabel('Imag. part of eigenvalues');
  
  exact_lambda
  computed_lambda
  legend('Exact eigenvalues','Computed eigenvalues')
  
  title(str_xlabel,'fontsize',12)
  
end
