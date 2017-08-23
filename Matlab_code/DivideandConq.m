% ----------------------------------------
% Computes algorithm of Divide-and-Conquer:
% eigenvalues will be roots of the secular equation and will lie
% on the diagonal of the output matrix L.
% In the output matrix Q will be corresponding eigenvectors.
% ----------------------------------------

function [Q,L] = DivideandConq(T)
  % Compute size of input matrix T:
  [m,n] = size(T);
  
  % here we will divide the matrix
  m2 = floor(m/2);
  
  %if m=0 we shall return
  if m2 == 0 %1 by 1
    Q = 1; L = T;
    return;
    %else we perform recursive computations
  else
    [T,T1,T2,bm,v] = formT(T,m2);
    
    %recursive computations
    [Q1,L1] = DivideandConq(T1);
    [Q2,L2] = DivideandConq(T2);
    
    %pick out the last and first columns of the transposes:
    Q1T = Q1';
    Q2T = Q2';
    u = [Q1T(:,end); Q2T(:,1)];
    
    %Creating the D-matrix:
    D = zeros(n);
    D(1:m2,1:m2) = L1;
    D((m2+1):end,(m2+1):end) = L2;
    
    % The Q matrix (with Q1 and Q2 on the "diagonals")
    Q = zeros(n);
    Q(1:m2,1:m2) = Q1;
    Q((m2+1):end,(m2+1):end) = Q2;
    
    %Creating the matrix B, which determinant is the secular equation:
    % det B = f(\lambda)=0
    B = D+bm*u*u';
    
    % Compute eigenvalues as roots of the secular equation
    %  f(\lambda)=0  using Newton's method
    eigs = NewtonMethod(D,bm,u);
    Q3 = zeros(m,n);
    
    % compute eigenvectors for corresponding eigenvalues
    for i = 1:length(eigs)
      Q3(:,i) = (D-eigs(i)*eye(m))\u;
      Q3(:,i) = Q3(:,i)/norm(Q3(:,i));
    end
    
    %Compute  eigenvectors of the original input matrix T
    Q = Q*Q3;
    
    % Present eigenvalues  of the original matrix input T
    %(they will be on diagonal)
    L = zeros(m,n);
    L(1:(m+1):end) = eigs;
    
    return;
  end
  
end

% Compute T1, T2  constant bm  and the vector v
%from the input matrix A.

function [T,T1,T2,bm,v] = formT(A,m)
  
  T1 = A(1:m,1:m);
  T2 = A((m+1):end,(m+1):end);
  bm = A(m,m+1);
  
  T1(end) = T1(end)-bm;
  T2(1) = T2(1)-bm;
  
  v = zeros(size(A,1),1);
  v(m:m+1) = 1;
  
  T = zeros(size(A));
  T(1:m,1:m) = T1;
  T((m+1):end,(m+1):end) = T2;
  
end

% compute eigenvalues in the secular equation
% using the Newton's method

function eigs = NewtonMethod(D,p,u)
  [m,n] = size(D);
  
  %The initial guess in  the Newton's method
  % will be the numbers d_i
  startingPoints = sort(diag(D));
  
  %if p > 0 we have an eigenvalue on the right, else on the left
  if p >= 0
    startingPoints = [startingPoints; startingPoints(end)+10000];
  elseif p < 0
    startingPoints = [startingPoints(1)-10000; startingPoints];
  end
  
  eigs = zeros(m,1);
  
  % tolerance in Newton's method
  convCriteria = 1e-05;
  
  % step in the approximation of the derrivative
  % in Newton's method
  dx = 0.00001;
  
  %plot the secular equation
  X = linspace(-3,3,1000);
  for t = 1:1000
    y(t) =SecularEqEval(D,p,u,X(t),m,n);
  end
  plot(X,y, 'LineWidth',2)
  axis([-3 3 -5 5])
  legend('graph of the secular equation $f(\lambda)=0$')
  
  %Start  Newton's method
  for i = 1:m
    %the starting value of lambda
    currentVal = (startingPoints(i)+startingPoints(i+1) )/ 2;
    
    % this value is used inthe stoppimg criterion below
    currentVal2 = inf;
    %  computed secular equation for \lambda=currentVal
    fCurr = SecularEqEval(D,p,u,currentVal,m,n);
    
    rands = 0;
    k =0;
    j = 0;
    
    if  ~((startingPoints(i+1)-startingPoints(i)) < 0.0001)
      while ~(abs(fCurr) < convCriteria)
        
        %compute value of the function  dfApprox with small step dx to
        %approximate derivative
        fval2 = SecularEqEval(D,p,u,currentVal+dx,m,n);
        fval1 = SecularEqEval(D,p,u,currentVal,m,n);
        dfApprox = (fval2-fval1)/dx;
        
        % compute new value of  currentVal in Newton's method,
        % or perform one iteration in Newton's method
        currentVal = currentVal - fCurr/dfApprox;
        
        % check: if we are  outside of the current range, reinput inside:
        if currentVal <= startingPoints(i)
          currentVal= startingPoints(i)+0.0001;
          k=k+1;
        elseif currentVal >= startingPoints(i+1);
          currentVal= startingPoints(i+1)-0.0001;
          k=k+1;
        elseif dfApprox == Inf || dfApprox == -Inf
          currentVal= startingPoints(i) + ...
          rand*(startingPoints(i+1)-startingPoints(i));
          rands = rands+1;
        end
        
        j=j+1;
        
        fCurr = SecularEqEval(D,p,u,currentVal,m,n);
        
        if k > 10 || j > 50;
          tempVec = [startingPoints(i),startingPoints(i+1)];
          [val,ind] = min(abs([startingPoints(i),startingPoints(i+1)]-currentVal));
          if ind == 1
            currentVal = tempVec(ind)+0.00001;
          else
            currentVal = tempVec(ind)-0.00001;
          end
          break;
        elseif currentVal2 == currentVal || rands > 5 || isnan(currentVal) || isnan(fCurr)
          currentVal = currentVal2;
          break;
        end
        %save last value:
        currentVal2 = currentVal;
      end
    end
    
    %assigning eigenvalue in the right order
    eigs(i) = currentVal;
    
  end
  
end

% evaluate the  secular equation in Newton's method for the computed
% eigenvalue x
function fVal = SecularEqEval(D,p,u,x,m,n)
  
  fVal = 1+p*u'*inv((D-x*eye(m,n)))*u;
  
end

