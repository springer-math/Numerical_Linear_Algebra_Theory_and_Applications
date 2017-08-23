% Polynomial evaluation using Horner's rule.
% Bisection algorithm and computation of error bounds.

clear
close all
clc

N = 8e3;
%x = linspace(1.85,2.15,N);
%x = linspace(1.87,2.13,N);
%x = linspace(1.92,2.08,N);
%x = linspace(1.85,2.15,N);
%x = linspace(-2.0,6.0,N);
%x = linspace(1.9,2.1,N);
%x = linspace(1.5,2.5,N);
%exact coefficients of polynomial p(x) = (x - 2)^9
%a = [-512 2304 -4608 5376 -4032 2016 -672 144 -18 1];

%exact coefficients of polynomial p(x) = (x - 9)^9
a = [-387420489 387420489 -172186884 44641044 -7440174 ...
826686 -61236 2916 -81 1];

% mesh for polynomial p(x)=(x-9)^9

%input interval for figures in the book

%x = linspace(-1.0,18,N);
%x = linspace(-3.0,15,N);
%x = linspace(8.5,9.5,N);
x = linspace(8.7,9.3,N);

%eps = 2.2204460492503131e-16; %machine epsilon in MATLAB
eps =  0.5e-16;

y_horner = zeros(N,1);

for i = 1:N
  
  [P,bp] = evaluate_polynomial_by_Horners_rule(a,x(i),eps);
  y_horner(i) = P;
  y_horner_upper(i) = P + bp;
  y_horner_lower(i) = P - bp;
  
end

%computation of exact polynomial
%y = (x - 2).^9;
y = (x - 9).^9;

%%%% run classical bisection: whithout computed error bound break   %%%%%

%  determine bounds for polynomial  (x - 9).^9
x_left = 8;  %left bound p(x_left) < 0
x_right = 11; %right bound p(x_right) > 0

[P,bp] = evaluate_polynomial_by_Horners_rule(a,x_left,eps);
p_left = P;
[P,bp] = evaluate_polynomial_by_Horners_rule(a,x_right,eps);
p_right = P;

%check

if p_left > 0 || p_right < 0
  disp('choose another interval')
  exit
end

iterations = 0;

while x_right - x_left > 2*eps
  
  iterations = iterations + 1;
  
  x_mid = (x_right + x_left)/2;
  [P,bp] = evaluate_polynomial_by_Horners_rule(a,x_mid,eps);
  p_mid = P;
  
  if p_left*p_mid < 0
    x_right = x_mid;
    p_right = p_mid;
  end
  if p_right*p_mid < 0
    x_left = x_mid;
    p_left = p_mid;
  end
  if p_mid == 0
    x_left = x_mid;
    x_right = x_mid;
  end
  
end

Root = (x_left + x_right)/2
% iterations

y_horner = zeros(N,1);

for i = 1:N
  
  [P,bp] = evaluate_polynomial_by_Horners_rule(a,x(i),eps);
  y_horner(i) = P;
  y_horner_upper(i) = P + bp;
  y_horner_lower(i) = P - bp;
  
  error(i) = abs(bp/P);
  
  log_error(i) = -log(abs(bp/P));
  
  % here we compute error between  computed and exact values of polynomial at x(i)
  ComputedErrors(i) = P- ((x(i)- 9).^9);
  
  LogCompEr(i) = -log(abs(ComputedErrors(i)./P));
  
end

figure(1)
plot(x,y_horner,'k.')
hold on
plot(x,y,'r','linewidth',2)
legend('Horners rule (8000 points)', 'exact p(x)')
xlabel('Input interval for x')

hold off

figure(2)
plot(x,log_error);

hold on
plot(x, LogCompEr, '. r');

legend('log error = -log(abs(bp/P))', '  -log(abs((P - p(x))/P)) ')

xlabel('input interval for x')

hold off

figure(3)

plot(x,error);

legend('error = abs(bp/P)')

xlabel('input interval for x')

figure(4)
plot(x,y,'k',x,y_horner_upper,'r--',x,y_horner_lower,'b--')
xlabel('input interval for x')
legend('exact  polynomial ','upper bound','lower bound')

