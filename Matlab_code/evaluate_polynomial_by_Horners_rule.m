function [P,bp] = evaluate_polynomial_by_Horners_rule(a,x,eps)
  % Parameters: a contains the coeficients of the plynomial  p(x)
  % P is the value of p(x) at x.
  % eps is the mechine epsilon
  d = numel(a);
  
  P = a(d);
  bp = abs(a(d));
  
  for i = d - 1:(-1):1
    
    P = x*P + a(i);
    bp = abs(x)*bp + abs(a(i));
    
  end
  
  %error bound
  bp = 2*(d - 1)*eps*bp;
  
end
