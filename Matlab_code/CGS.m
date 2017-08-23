% ----------------------------------------
%            Classical Gram-Schmidt (CGS) orthogonalization process
%  and solution of the linear least square problem  using CGS.
% ----------------------------------------

% size of our matrix A is m-by-n
m= 6;
n=3;

% vector of the right hand side
y=zeros(1,m);

A=[1,0,0;
0,1,0;
0,0,1;
-1, 1,0;
-1,0,1;
0,-1,1];

y = [1237,1941,2417,711,1177,475];

% allocate matrices q and r for QR decomposition

q=[];
r=[];

%QR decomposition using classical Gram-Schmidt orthogonalization
for k=1:1:n
  q(:,k)=A(:,k);
  for j=1:1:k-1
    r(j,k)=q(:,j)'*A(:,k);
    q(:,k)=q(:,k)-r(j,k)*q(:,j);
  end
  r(k,k)=norm(q(:,k));
  q(:,k)=q(:,k)/r(k,k);
end

%compute solution of the system Ax = QR x =  y
% by backward substitution:  R x = Q^T y

b=[];

% compute right hand side Q^T y
b=q'*y';

% perform backward substitution to get solution x = R^(-1)  Q^T y
% obtain solution in b
for i=n:-1:1
  for k=n:-1:i+1
    b(i)=b(i)-b(k)*r(i,k);
  end
  b(i)=b(i)/r(i,i);
end

