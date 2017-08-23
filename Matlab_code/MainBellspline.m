% ----------------------------------------
%   Solution of least squares  problem  min_x || Ax - y ||_2
%   using the method of normal equations, QR decomposition
%   and SVD decomposition.
%   Matrix  A is constructed using bellsplines.
%   Program performs fitting to the function y = sin(pi*x/5) + x/5
% ----------------------------------------

clc
clear
clf
close all
format short

% input interval on which we fit the function
interval=10;

% junction points

T=linspace(-10,interval,7)';

% Define number of measurement points m
m=30;
x=linspace(-10,interval,m)';

%exact function to be fitted
b=sin(pi*x/5) +x/5;

% construct matrix A with bellsplines
%Number of bellsplines should be number of junction points +2

A=fbell(x,T);

%solution of system Ax = b using different methods for solution
% of least squares problem.
tic
% use method of normal equations
xHatChol = LLSChol(A,b);
toc
tic
%use SVD decomposition of A
xHatSVD = LLSSVD(A,b);
toc
tic
% use QR decomposition of A
xHatQR = LLSQR(A,b);
toc

% compute condition number of A
cond(A)

% use iterative refinement of the obtained solution
%  via Newton's method
% choose tolerance in Newton's method
tol =0.2;

y=  newtonIR(A,xHatChol,b,tol);
y1= newtonIR(A,xHatQR,b,tol);
y2= newtonIR(A,xHatSVD,b,tol);

% compute relative errors

eC=norm(A*xHatChol-b)/norm(b);
eS=norm(A*xHatSVD-b)/norm(b);
eQ=norm(A*xHatQR-b)/norm(b);

disp(' --------------Computed relative errors ------------------- ')
disp('      Method of normal eq.           QR              SVD')
disp('')

disp([eC eS eQ ])

disp('Computed relative errors after iterative refinement via Newton method ')
disp('    Method of normal eq.            QR             SVD')
disp('')

disp([norm(A*y-b)/norm(b) norm(A*y1-b)/norm(b) norm(A*y2-b)/norm(b)])

% Plot results

figure(1)
%plot(t,A,'linewidth',2)

plot(x,A,'linewidth',2)

m =size(A,2);
str_xlabel = [' number of bellsplines=', num2str(m)];
title(str_xlabel)

figure('Name','Cholesky')
title('Cholesky')
plot(x,b,'o- r', 'linewidth',2)
hold on
plot(x,A*xHatChol,' *- b', 'linewidth',2)
legend('exact  ', 'B-spline degree 3, Cholesky');

figure('Name','QR')
plot(x,b,'o- r', 'linewidth',2)
hold on
plot(x,A*xHatQR,'* - b', 'linewidth',2)
legend('exact ', 'B-spline degree 3, QR');

figure('Name','SVD')
title('SVD')
plot(x,b,'o- r', 'linewidth',2)
hold on
plot(x,A*xHatSVD,'*- b', 'linewidth',2)
legend('exact ', 'B-spline degree 3, SVD');
