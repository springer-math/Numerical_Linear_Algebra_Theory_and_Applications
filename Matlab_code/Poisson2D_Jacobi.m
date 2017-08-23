% ----------------------------------------
% Main program for the solution of Poisson's equation
%  - a laplace = f in 2D using iterative Jacobi method
% ----------------------------------------

close all
clc
clear
clf
%Define input parameters
n=20; % number of inner nodes in one direction.
a_amp = 12; %  amplitude for the function a(x_1,x_2)
f_amp = 1; %  we can choose f=1, 50, 100
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

% If f is constant.
% f = f_amp*ones(n^2,1);

% If f is Gaussian function.
f=zeros(n^2,1);
for i=1:n
  for j=1:n
    f(n*(i-1)+j)= f_amp*exp(-((i*h-x_0)^2/(2*c_x^2)...
    +(j*h-y_0)^2/(2*c_y^2)));
  end
end

%  Compute vector of right hand side
%  b = D^(-1)*f   computed as b(i,j)=f(i,j)/a(i,j)

b=zeros(n^2,1);
for i=1:n
  for j=1:n
    b(n*(i-1)+j)= f(n*(i-1)+j)/C(i,j); % Use coefficient matrix C or
    % diagonal matrix D to get a(i,j)
  end
end

% ----------------------------------------
% ---  Solution of 1/h^2 S*u = b using Jacobi's method, version I
% ----------------------------------------

err = 1;  k=0;  tol=10^(-9);

w_old = ones(length(S),1);
L=tril(S,-1);
U=L';
Dinv=diag(diag(S).^(-1));
R=Dinv*(-L-U);
c=Dinv*h^2*b;

while(err>tol)
  w_new = R*w_old +c;
  k=k+1;
  
  % stopping criterion: choose one of two
  err = norm(w_new-w_old);
  %  err = norm(S*w_new - h^2*b);
  w_old = w_new;
end

disp('-- Number of iterations in the version I of Jacobi method ----------')
k

% ----------------------------------------
% Plots and figures for version I
% ----------------------------------------

% sort the data in u into the mesh-grid, the boundary nodes are zero.
V_new = zeros(n+2,n+2);
for i=1:n
  for j=1:n
    V_new(i+1,j+1) = w_new(j+n*(i-1));
  end
end

%% plotting
x1=0:h:1;
y1=0:h:1;

figure(1)

subplot (2,2,1)
surf(x1,y1,V_new) % same plot as above, (x1, y1 are vectors)
view(2)
colorbar
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')
title( ['solution u(x_1,x_2) Jacobi version I ',...
',  N = ',num2str(n),', iter. = ',num2str(k)])

subplot (2,2,2)

surf(x1,y1,V_new) % same plot as above
colorbar
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')

title( ['solution u(x_1,x_2) Jacobi version I',...
',  N = ',num2str(n),', iter. = ',num2str(k)])

% Plotting a(x,y)
Z_a= zeros(n+2);
for i=1:(n+2)
  for j=1:(n+2)
    Z_a(i,j)= 1 + a_amp*exp(-((i*h-x_0)^2/(2*c_x^2)...
    +(j*h-y_0)^2/(2*c_y^2)));
  end
end

subplot (2,2,3)
surf(x1,y1,Z_a)
xlabel('x_1')
ylabel('x_2')
zlabel('a(x_1,x_2)')
title( ['coefficient a(x_1,x_2) with A = ',num2str(a_amp)])

%  plott the function f(x,y)
Z_f= zeros(n+2);
for i=1:(n+2)
  for j=1:(n+2)
    Z_f(i,j)=f_amp*exp(-((x1(i)-x_0)^2/(2*c_x^2)...
    +(y1(j)-y_0)^2/(2*c_y^2)));
  end
end

subplot (2,2,4)
surf(x1,y1,Z_f)
xlabel('x_1')
ylabel('x_2')
zlabel('f(x_1,x_2)')
title( ['f(x_1,x_2) with A_f = ',num2str(f_amp)])

% ----------------------------------------
% ---  Jacobi's method, version II --------------------------
% ----------------------------------------

k=0; err = 1;
V_old = zeros(n,n);
V_new = zeros(n,n);
F=vec2mat(b,n)';
X=diag(ones(1,n-1),-1);
X=X+X';

while(err>tol)
  V_new = (X*V_old + V_old*X' + h^2*F)/4;
  k=k+1;
  err = norm(V_new-V_old);
  V_old = V_new;
end

%apply boundary conditions
V_new = [zeros(1,n+2); zeros(n,1) V_new zeros(n,1);zeros(1,n+2)]

disp('-- Number of iterations in the version II of Jacobi method ----------')
k

figure(2)
% ----------------------------------------
% Plots and figures for version II
% ----------------------------------------

%% plotting
x1=0:h:1;
y1=0:h:1;

subplot (1,2,1)
surf(x1,y1,V_new) % same plot as above, (x1, y1 are vectors)
view(2)
colorbar
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')
title( ['solution u(x_1,x_2) Jacobi version II ',...
',  N = ',num2str(n),', iter. = ',num2str(k)])

subplot (1,2,2)

surf(x1,y1,V_new) % same plot as above
colorbar
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')

title( ['solution u(x_1,x_2)  Jacobi version II',...
',  N = ',num2str(n),', iter. = ',num2str(k)])

% ----------------------------------------
% ---  Jacobi's method, version III --------------------------
% ----------------------------------------

err = 1;  k=0;  tol=10^(-9);
% Initial guess
uold = zeros(n+2, n+2);
unew= uold;

% counter for iterations
k = 0;

while(err > tol)
  for i = 2:n+1
    for j = 2:n+1
      unew(i, j) = (uold(i-1, j) + uold(i+1, j) + uold(i, j-1) + uold(i, j+1)  + h^2*b(n*(i-2)+j-1))/4.0;
    end
  end
  k = k+1;
  err = norm(unew-uold);
  uold = unew;
end

u = reshape(unew(2:end-1, 2:end-1)', n*n, 1);

disp('-- Number of iterations in the version III of Jacobi method ----------')

k

figure(3)
% ----------------------------------------
% Plots and figures for version III
% ----------------------------------------

% sort the data in u into the mesh-grid, the boundary nodes are zero.
V_new = zeros(n+2,n+2);
for i=1:n
  for j=1:n
    V_new(i+1,j+1) = u(j+n*(i-1));
  end
end

%% plotting
x1=0:h:1;
y1=0:h:1;

subplot (1,2,1)
surf(x1,y1,V_new) % same plot as above, (x1, y1 are vectors)
view(2)
colorbar
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')
title( ['solution u(x_1,x_2) Jacobi version III ',...
',  N = ',num2str(n),', iter. = ',num2str(k)])

subplot (1,2,2)

surf(x1,y1,V_new) % same plot as above
colorbar
xlabel('x_1')
ylabel('x_2')
zlabel('u(x_1,x_2)')

title( ['solution u(x_1,x_2)  Jacobi version III',...
',  N = ',num2str(n),', iter. = ',num2str(k)])

