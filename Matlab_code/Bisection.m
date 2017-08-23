% ----------------------------------------
% Find all eigenvalues of the matrix A ion the input interval [a,b)
% ----------------------------------------

% define size n of the  n-by-n matrix A
n=5;

% Generate the symmetric tridiagonal  matrix A
A=randomTridiag(n);

% Set bounds for the interval [a,b) in the algorithm and the tolerance
a=-100;b=100;
tol=0.000001;

%Define functions for the worklist

DeleteRowInWorklist=@(Worklist,linenr) ChangeRowInWorklist(Worklist,linenr,'delete');
InsertRowInWorklist=@(Worklist,LineToAdd)...
ChangeRowInWorklist(Worklist,LineToAdd,'add');

% Set the info for the first worklist
na=Negcount(A,a);
nb=Negcount(A,b);
Worklist=[];

%If no eigenvalues are found on the interval [a,b) then save an empty worklist
if na ~= nb
  Worklist=InsertRowInWorklist(Worklist,[a,na,b,nb]);
end

while numel(Worklist) ~= 0
  [Worklist, LineToWorkWith ]= DeleteRowInWorklist(Worklist,1);
  
  low=LineToWorkWith(1);
  n_low=LineToWorkWith(2);
  up=LineToWorkWith(3);
  n_up=LineToWorkWith(4);
  
  % if the upper and lower bounds are close enough we  print out this interval
  if (up-low)< tol
    NrOfEigVal = n_up-n_low;
    fprintf('We have computed %3.0f eigenvalues in the interval [%4.4f,%4.4f) \n', ...
    NrOfEigVal,low,up);
  else
    % Perform the bisection step
    mid= (low+up)/2;
    n_mid= Negcount(A,mid);
    if n_mid > n_low
      Worklist = InsertRowInWorklist(Worklist,[low,n_low,mid,n_mid]);
    end
    if n_up>n_mid
      Worklist = InsertRowInWorklist(Worklist,[mid,n_mid,up,n_up]);
    end
  end
end

