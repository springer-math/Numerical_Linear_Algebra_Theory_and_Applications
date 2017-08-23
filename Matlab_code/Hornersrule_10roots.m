% Polynomial p(x) = (x -1)^2(x-2)(x-3)(x-4)(x-5)(x-7)(x-9)(x-11)(x-15)(x-17)
% Polynomial evaluation using Horner's rule.
% Computation of error bounds.

clear
close all
clc

N = 8e3;

%  exact coefficients of polynomial
%  p(x) = (x -1)^2(x-2)(x-3)(x-4)(x-5)(x-7)(x-9)(x-11)(x-15)(x-17)

%input interval for  p(x)

x = linspace(-1.0, 20,N);

%Get coefficients of   polynomial p(x) = (x -1)^2(x-2)(x-3)(x-4)(x-5)
syms t;

a=double(coeffs((t-1).^2*(t-2)*(t-3)*(t-4)*(t-5)*(t-7)*(t-9)*(t-11)*(t-15)*(t-17),t))

%eps = 2.2204460492503131e-16; %machine epsilon in MATLAB
eps =  0.5e-16;

y_horner = zeros(N,1);

for i = 1:N
  
  [P,bp] = evaluate_polynomial_by_Horners_rule(a,x(i),eps);
  y_horner(i) = P;
  y_horner_upper(i) = P + bp;
  y_horner_lower(i) = P - bp;
  
end

y_horner = zeros(N,1);

y=0;

for i = 1:N
  
  [P,bp] = evaluate_polynomial_by_Horners_rule(a,x(i),eps);
  y_horner(i) = P;
  y_horner_upper(i) = P + bp;
  y_horner_lower(i) = P - bp;
  
  error(i) = abs(bp/P);
  
  log_error(i) = -log(abs(bp/P));
  
  % here we compute error between  computed and exact values of polynomial at x(i)
  ComputedErrors(i) = P- (x(i)-1).^2*(x(i)-2)*(x(i)-3)*(x(i)-4)*(x(i)-5)*(x(i)- 7)*
  (x(i) - 9)*(x(i) - 11)*(x(i) -15)*(x(i)-17);
  
  y(i) =  (x(i)-1).^2*(x(i)-2)*(x(i)-3)*(x(i)-4)*(x(i)-5)*(x(i)- 7)*(x(i) - 9)*(x(i) -11)*
  (x(i) -15)*(x(i)-17);
  
  LogCompEr(i) = -log(abs(ComputedErrors(i)./P));
  
end

figure(1)
plot(x,y_horner,'k.')
hold on
plot(x,y,'r','linewidth',2)
legend('Horners rule (8000 points)', 'exact p(x) ')
xlabel('Input interval for x')

hold off

figure(2)
plot(x,log_error);

hold on
plot(x, LogCompEr, '. r');

legend('  estimated bound', '  computed bound ')

xlabel('input interval for x')

hold off

figure(3)
plot(x,y,'k',x,y_horner_upper,'r--',x,y_horner_lower,'b--')
xlabel('input interval for x')
legend('exact p(x)','upper bound','lower bound')

