% ----------------------------------------
%   Solution of least squares  problem  min_x || Ax - y ||_2
%   using the method of normal equations.
%   Matrix  A is constructed as a Vandermonde matrix.
%   Program performs fitting to the function y = sin(pi*x/5) + x/5
% ----------------------------------------

d=5;  % degree of the polynomial
m=10;%number of discretization points or rows in the matrix A

x=zeros(1,m);
y=zeros(1,m);
A=[];
for i=1:1:m
  x = linspace(-10.0,10.0,m);
  %  exact function which should be approximated
  y(i)= sin(pi*x(i)/5) + x(i)/5;
end

% construction of a Vamdermonde matrix

for i=1:1:m
  for j=1:1:d+1
    A(i,j)=power(x(i),j-1);
  end
end

% computing the right hand side in the method of normal equations
c=A'*y';

% computing matrix in the left hand side in the method of normal equations
C=A'*A;

l=zeros(d+1);

% solution of the normal equation using Cholesky decomposition

for j=1:1:d+1
  s1=0;
  for k=1:1:j-1
    s1=s1+l(j,k)*l(j,k);
  end
  l(j,j)=(C(j,j)-s1)^(1/2);
  for i=j+1:1:d+1
    s2=0;
    for k=1:1:j-1
      s2=s2+l(i,k)*l(j,k);
    end
    l(i,j)=(C(i,j)-s2)/l(j,j);
  end
end
for i=1:1:d+1
  for k=1:1:i-1
    c(i)=c(i)-c(k)*l(i,k);
  end
  c(i)=c(i)/l(i,i);
end
for i=d+1:-1:1
  for k=d+1:-1:i+1
    c(i)=c(i)-c(k)*l(k,i);
  end
  c(i)=c(i)/l(i,i);
end

figure(1)
plot(x,y,'o- r', 'linewidth',1)
hold on

% compute approximation to this exact polynomial with comp. coefficients c

approx = A*c;
plot(x,approx,'*- b', 'linewidth',1)
hold off

str_xlabel = ['poly.degree d=', num2str(d)];

legend('exact  sin(pi*x(i)/5) + x(i)/5',str_xlabel);

xlabel('x')

% computation of the relative error as
%  norm(approx. value - true value) / norm(true value)
e1=norm(y'- approx)/norm(y')
