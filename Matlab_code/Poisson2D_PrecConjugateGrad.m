% ----------------------------------------
% Main program for the solution of Poisson's equation
%  - a laplace = f in 2D using  Preconditioned Conjugate Gradient Method
% ----------------------------------------

close all
%Define input parameters
n=20; % number of inner nodes in one direction.
a_amp = 12; %  amplitude for the function a(x_1,x_2)
f_amp = 1; % we can set f = 1, 50, 100
x_0=0.5;
y_0=0.5;
c_x=1;
c_y=1;

h = 1/(n+1); % define step length

% ----------------------------------------
% Computing all matrices and vectors
% ----------------------------------------
% Generate a n*n by n*n stiffness matrix
S = DiscretePoisson2D(n);

%% generate coefficient matrix of a((x_1)_i,(x_2)_j) = a(i*h,j*h)
C = zeros(n,n);
for i=1:n
  for j=1:n
    C(i,j) = 1 + a_amp*exp(-((i*h-x_0)^2/(2*c_x^2)...
    +(j*h-y_0)^2/(2*c_y^2)));
  end
end
% create diagonal matrix from C
D = zeros(n^2,n^2);
for i=1:n
  for j=1:n
    D(j+n*(i-1),j+n*(i-1)) = C(i,j);
  end
end

%% calculate load vector f

% If f is constant.
% f = f_amp*ones(n^2,1);

% If f is Gaussian function.
f=zeros(n^2,1);
for i=1:n
  for j=1:n
    f(n*(i-1)+j)=f_amp*exp(-((i*h-x_0)^2/(2*c_x^2)...
    +(j*h-y_0)^2/(2*c_y^2)));
  end
end

% 1. Compute vector of right hand side
%  b = D^(-1)*f given by b(i,j)=f(i,j)/a(i,j)

b=zeros(n^2,1);
for i=1:n
  for j=1:n
    b(n*(i-1)+j)=f(n*(i-1)+j)/C(i,j); % Use coefficient matrix C or
    % diagonal matrix D to get a(i,j)
  end
end

% ----------------------------------------
% --- Preconditioned conjugate gradient method (PCGM):
%  choose  different preconditioners:
% Cholesky factorization, Jacobi  preconditioner, block Jacobi preconditioner
% ----------------------------------------
% We now have system to solve: 1/h^2 S u = b

%initialize preconditioner
Ssparse = sparse(S);

% Preconditioner: preconditioner matrix  here is incomplete
% Cholesky factorization of S

cond = ichol(Ssparse); cond=cond*cond'; cond=full(inv(cond));

% Preconditioner: preconditioner matrix  here is
% Jacobi preconditioner.
% Results are the same as in usual conjugate gradient update

%M = diag(diag(S));
%cond = diag(1.0./diag(M));

% Preconditioner: preconditioner matrix  here is
%Block Jacobi Preconditioner
%blockSize = 2; % size of blocks
%cond = zeros(n^2);
%Iinds = ceil( (1:(blockSize*n^2))/blockSize);
%Jinds = blockSize*ceil( (1:(blockSize*n^2))/blockSize^2)-(blockSize-1) ...
%       + repmat%(0:blockSize-1,1,n^2);
%vecInds = sub2ind(size(S),Iinds, Jinds);
%cond(vecInds) = S(vecInds);

%initialize parameters in the method
err = 1; x=0; r0= h^2*b; p=cond*h^2*b; y0=cond*r0; tol=10^(-9);
k=0;

while(err>tol)
  z = S*p;
  nu = (y0'*r0)/(p'*z);
  x = x + nu*p;
  r1 = r0 - nu*z;
  y1 = cond*r1;
  mu = (y1'*r1)/(y0'*r0);
  p = y1 + mu*p;
  r0=r1;
  y0=y1;
  err = norm(r0);
  k=k+1;
end

disp('-- Number of iterations in  Preconditioned conjugate gradient method (PCGM) ')
k

% ----------------------------------------
% Plots and figures.
% ----------------------------------------

% sort the data in u into the mesh-grid, the boundary nodes are zero.
Z = zeros(n+2,n+2);
for i=1:n
  for j=1:n
    Z(i+1,j+1) = x(j+n*(i-1));
  end
end

%% plotting
x1=0:h:1;
y1=0:h:1;

subplot(2,2,1)
surf(x1,y1,Z) % same plot as above, (x1, y1 are vectors)
view(2)
colorbar
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')
title( ['u(x_1,x_2) in PCGM',...
',  N = ',num2str(n)])

subplot(2,2,2)
surf(x1,y1,Z) % same plot as above
colorbar
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')
title( ['u(x_1,x_2) in PCGM ',...
', N = ',num2str(n)])

% Plotting a(x,y)
Z_a= zeros(n+2);
for i=1:(n+2)
  for j=1:(n+2)
    Z_a(i,j)= 1 + a_amp*exp(-((i*h-x_0)^2/(2*c_x^2)...
    +(j*h-y_0)^2/(2*c_y^2)));
  end
end

subplot(2,2,3)
surf(x1,y1,Z_a)
xlabel('x_1')
ylabel('x_2')
zlabel('a(x_1,x_2)')
title( ['a(x_1,x_2) with A = ',num2str(a_amp)])

%  plott the function f(x,y)
Z_f= zeros(n+2);
for i=1:(n+2)
  for j=1:(n+2)
    Z_f(i,j)=f_amp*exp(-((x1(i)-x_0)^2/(2*c_x^2)...
    +(y1(j)-y_0)^2/(2*c_y^2)));
  end
end

subplot(2,2,4)
surf(x1,y1,Z_f)
xlabel('x_1')
ylabel('x_2')
zlabel('f(x_1,x_2)')
title( ['f(x_1,x_2) with A_f = ',num2str(f_amp)])

