% ----------------------------------------
%              Method of  QR iteration with  shift  sigma=A(n,n)
% ----------------------------------------
clc
%clear all
%close all
eps = 1e-09;
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
  %choose shift
  %sigma=1.0;
  sigma=A(n,n);
  %sigma=A(1,1);
  
  % get  exact eigenvalues in sorted order
  exact_lambda = sort(eig(A));
  %%% Method of QR iteration with shift
  
  for k = 1:100
    
    A = A - sigma*eye(n);
    
    [Q,R] = qr(A);
    % end
    
    A = R*Q + sigma*eye(n);
    
    %compute shift
    sigma=A(n,n);
    
    % %%%%%%%%% Find eigenvalues from Real Schur block
    j =2; count =1;
    eigs = zeros(1,n);
    while (j <=n)
      %real eigenvalues
      if(abs(A(j,j-1)) < 1e-7)
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
    %*******************************
    
    computed_lambda = sort(eigs)';
    
    if(norm(abs(computed_lambda - lambda0 ))<eps )
      break ;
    end
    lambda0 = computed_lambda ;
    iter = iter + 1;
  end
  %***************************************
  str =['Comp. eig.:' num2str(computed_lambda')];
  str=[str, ', Ex. eig.:' num2str(exact_lambda',2)];
  str_xlabel = ['Example ',num2str(i), '.  Nr.it. in  QR it. with shift:', num2str(iter)];
  
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
