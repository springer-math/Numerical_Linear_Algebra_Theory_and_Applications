% ----------------------------------------
%   Matrix  B is constructed using  bellsplines.
%   Input arguments: T - column vector with junction points,
%   x are measurement ponts (discretization points).
% ----------------------------------------

function B=fbell(x,T)
  
  m=length(x);
  N=length(T);
  epsi=1e-14;
  
  %construct  N+6 column vector
  a=[T(1)*[1 1 1]'; T; T(N)*(1+epsi)*[1 1 1]'];
  n=N+5;
  C=zeros(m,n);
  for k=1:n
    I=find(x>=a(k) & x<a(k+1));
    if ~isempty(I)
      C(I,k)=1;
    end
  end
  for j=1:3
    B=zeros(m, n-j);
    for k=1:n-j
      d1=(a(k+j)-a(k));
      if abs(d1)<=epsi
        d1=1;
      end
      d2=(a(k+j+1)-a(k+1));
      if abs(d2)<=epsi
        d2=1;
      end
      B(:,k)=(x-a(k)).*C(:,k)/d1 + (a(k+j+1)-x).*C(:,k+1)/d2;
    end
    C=B;
  end
  
