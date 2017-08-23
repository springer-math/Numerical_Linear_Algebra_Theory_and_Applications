% main program for the solution of Poisson's equation
% - a laplace = f in 2D using Cholesky decomposition

close all
%Define input parameters
n=20; % number of inner nodes in one direction.

A_1 = 10; %  amplitude 1 for the rhs
A_2 = 10; %  amplitude 2 for the rhs

h = 1/(n+1); % define step length

% ----------------------------------------
% Computing all matrices and vectors
% ----------------------------------------
% Generate a n*n by n*n stiffness matrix
S = DiscretePoisson2D(n);

% factorize A=L*L^T using Cholesky decomposition
[L]=Cholesky(S);

%% generate coefficient matrix of a((x_1)_i,(x_2)_j) = a(i*h,j*h)
C = zeros(n,n);
for j=1:n
  for i=1:n
    C(i,j) = 1;
  end
end

%% compute load vector f

f=zeros(n^2,1);
for j=1:n
  for i=1:n
    f(n*(i-1)+j)= A_1*exp(-((i*h-0.25)^2/0.02...
    +(j*h-0.25)^2/0.02))+ A_2*exp(-((i*h-0.75)^2/0.02...
    +(j*h-0.75)^2/0.02));
  end
end

% ----------------------------------------
% Solving the linear system of equations using Gaussian elimination
% ----------------------------------------
% We have system A u = 1/h^2 (C*L*L^T) u = f

% 1. Compute vector of right hand side
% as b(i,j)=f(i,j)/a(i,j)

b=zeros(n^2,1);
for j=1:n
  for i=1:n
    b(n*(i-1)+j)=f(n*(i-1)+j)/C(i,j); % Use coefficient matrix C
    
  end
end

% We now have system to solve: 1/h^2 A u = b
% Use first LU decomposition: 1/h^2 (L L^T)  u = b
% 2. Compute v = L^(-1)*b by forward substitution.

v=ForwSub(L,b);

% We now have system 1/h^2 L^T u = v
% 3. Compute w = L^T^(-1)*v by backward substitution.

w=BackSub(L',v);

% 4. We now have system 1/h^2 u = w
% Compute finally solution as:  u=h^2*w
u=h^2*w;

% ----------------------------------------
% Plots and figures.
% ----------------------------------------

% sort the data in u into the mesh-grid, the boundary nodes are zero.
Z = zeros(n+2,n+2);
for j=1:n
  for i=1:n
    Z(i+1,j+1) = u(n*(i-1)+j);
  end
end

%% plotting
x1=0:h:1;
y1=0:h:1;

figure(1)
surf(x1,y1,Z) % same plot as above, (x1, y1 are vectors)
view(2)
colorbar
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')
title( ['u(x_1,x_2) with N = ',num2str(n)])

figure(2)
surf(x1,y1,Z) % same plot as above
colorbar
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')
title( ['u(x_1,x_2) with N = ',num2str(n)])
