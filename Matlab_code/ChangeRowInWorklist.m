% ----------------------------------------
% Add or remove rows to the  WorkList
% If action = 'add' then add a line to the Worklist, return Worklist and
%               new line
%   If action = 'delete' then delete the given line from  the Worklist, return
%               Worklist and deleted line
% ----------------------------------------

function [ Worklist , LineInQuestion] = ChangeRowInWorklist(Worklist,LINE,action)
  
  if strcmp(action,'delete')
    if (length(Worklist(:,1)) == 1)
      LineInQuestion=Worklist;
      Worklist=[];
    elseif (LINE==length(Worklist(:,1)))
      LineInQuestion = Worklist(LINE,:);
      Worklist=Worklist(1:(end-1),:);
    elseif (LINE==1)
      LineInQuestion = Worklist(LINE,:);
      Worklist=Worklist(2:end,:);
    else
      LineInQuestion = Worklist(LINE,:);
      Worklist=[Worklist(1:(LINE-1),:);Worklist((LINE+1):end,:)];
    end
    
  elseif strcmp(action,'add')
    LineInQuestion = LINE;
    if (length(Worklist) == 0)
      Worklist=LINE;
    else
      Worklist = [Worklist;LINE];
    end
    
  else
    fprintf('The third argument must be either delete or add!')
    
  end
end

