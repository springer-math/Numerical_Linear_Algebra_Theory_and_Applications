% ----------------------------------------
%   Solution of least squares  problem  min_x || Ax - y ||_2
% using QR decomposition. QR decomposition is performed via classical
%  Gram-Schmidt  (CGM) orthogonalization procedure.
%   Matrix  A is constructed as a Vandermonde matrix.
%   Program performs fitting to the function y = sin(pi*x/5) + x/5
% ----------------------------------------

d=5;  % degree of polynomial
m=10;%number of discretization points or rows in the matrix A
p=ones(1,d+1);
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

q=[];
r=[];

%QR decomposition via CGM

for i=1:1:d+1
  q(:,i)=A(:,i);
  for j=1:1:i-1
    r(j,i)=q(:,j)'*A(:,i);
    q(:,i)=q(:,i)-r(j,i)*q(:,j);
  end
  r(i,i)=norm(q(:,i));
  q(:,i)=q(:,i)/r(i,i);
end
b=[];
b=q'*y';
for i=d+1:-1:1
  for k=d+1:-1:i+1
    b(i)=b(i)-b(k)*r(i,k);
  end
  b(i)=b(i)/r(i,i);
end

figure(1)
plot(x,y,'o- r', 'linewidth',1)
hold on

% compute approximation to this exact polynomial with comp. coefficients b

approx = A*b;
plot(x,approx,'*- b', 'linewidth',1)
hold off

str_xlabel = ['poly.degree d=', num2str(d)];

legend('exact  sin(pi*x(i)/5) + x(i)/5',str_xlabel );

xlabel('x')

% computation of the relative error as
%  norm(approx. value - true value) / norm(true value)
e1=norm(y'- approx)/norm(y')
