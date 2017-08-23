% ----------------------------------------
%   Solution of the least squares problem $ \min_x || Ax - y ||_2 $
%   using the method of normal equations, QR decomposition
%   and SVD decomposition.
%   Matrix $A$ is constructed using linear splines.
%   The program performs fitting to the function $y = \sin(\pi*x/5) + x/5 $
% ----------------------------------------

clc
clear
clf
format long
close all

% Define the number of measurements or data points.
% It is also the number of columns in matrix $A$.
m = 100;

% the number of junction points
K = 5;

x = linspace(-10, 10.0, m)';
T = linspace(-10, 10.0, K)';

% function which we want to fit
b = sin(pi*x/5) + x/5;

A = zeros(m, K);

% construct matrix A using linear splines
for k = 1:K
  A(:,k) = fihatt(k, x, T);
end
% compute condition number of A
cond(A)

% solution of linear system $Ax = b$ by different methods

% using method of normal equations
xHatChol = LLSChol(A, b);

% using QR decomposition of A
xHatQR = LLSQR(A, b);

% using SVD decomposition of A
xHatSVD = LLSSVD(A, b);

disp(' Computed relative error ')
disp('      Method of normal eq.           QR              SVD')
disp('')

disp([norm(A*xHatChol-b)/norm(b) norm(A*xHatQR-b)/norm(b)   ...
norm(A*xHatSVD-b)/norm(b)])

% Method of iterative refinement via Newton's method

tol = 0.07;
refinedC = newtonIR(A, xHatChol, b, tol);
refinedQ = newtonIR(A, xHatQR, b, tol);
refinedS = newtonIR(A, xHatSVD, b, tol);

disp('Computed relative error after iterative refinement via Newton method ')
disp('    Method of normal eq.            QR             SVD')
disp('')

disp([norm(A*refinedC-b)/norm(b) norm(A*refinedQ-b)/norm(b) norm(A*refinedS-b)/norm(b)])

% Plot exact and computed functions

% choose the number of points to plot solution
x = linspace(-10, 10.0, 100)';
b = (sin(pi*x/5) + x/5);
A = zeros(100, K);

for k = 1:K
  A(:,k) = fihatt(k, x, T);
end

% choose method to be plotted

%method = 'cholesky';
%method = 'refinedcholesky';
%method = 'qr';
%method = 'refinedqr';
%method = 'svd';
method = 'refinedsvd';

switch lower(method)
          case 'cholesky'
  % Here, A is constructed by linear splines, approximated function is computed
  % via  the method of normal equations (Cholesky decomposition)
  solution = A*xHatChol;
case 'refinedcholesky'
  % Here, A is constructed by linear splines, approximated function is computed
  % via iterative refinement of the Cholesky-solution through the Newton method
  solution = A*refinedC;
case 'qr'
  % Here, A is constructed by linear splines, approximated function is computed
  % via QR decomposition
  solution = A*xHatQR;
case 'refinedqr'
  % Here, A is constructed by linear splines, approximated function is computed
  % via iterative refinement of the QR-solution through the Newton method
  solution = A*refinedQ;
case 'svd'
  % Here, A is constructed by linear splines, approximated function is computed
  % via SVD decomposition
  solution = A*xHatSVD;
case 'refinedsvd'
  % Here, A is constructed by linear splines, approximated function is computed
  % via iterative refinement of the SVD-solution through the Newton method
  solution = A*refinedS;
otherwise
  disp('Unknown method')
end

figure (1)
plot(x, b, 'o r', 'linewidth', 1)
hold on
plot(x, solution, ' - b', 'linewidth', 2)
legend('function', 'approx');
figure('Name', 'Hat functions')
plot(x, A, 'k')
