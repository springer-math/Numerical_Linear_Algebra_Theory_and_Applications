
% ----------------------------------------
% Run Hager's algorithm.
% ----------------------------------------

function [LowerBound] = HagersAlg(B)
  
  x=(1/length(B))*ones(length(B),1);
  
  iter=1;
  while iter < 1000
    w=B*x; xi=sign(w); z = B'*xi;
    if max(abs(z)) <= z'*x
      break
    else
      x= (max(abs(z))== abs(z));
    end
    iter = iter + 1;
  end
  LowerBound = norm(w,1);
end

