% ----------------------------------------
%  Iterative refinement using Newton's method.
%   Matrix A is m-by-n, m > n, the vector of the rhs b is of the size n.
% ----------------------------------------

function w=newtonIR(A,x,b,tol)
  
  relative_error=1;
  iter = 0;
  
  while relative_error > tol
    
    %compute residual
    r = A*x-b;
    d=A\r;
    x=x-d;
    iter = iter+1
    relative_error = norm(A*x - b)/norm(b)
    
    % here we introduce the maximal number of iterations
    % in Newton's method: if the relative error
    % is not rediced -  we terminate computations
    
    if iter  > 100
      break
    end
  end
  w=x;
