% Finds some feasable intercept ito. equation indeces, and it's
% coordinates, from equations e{}, conditions{}, start search index 
% i = [eq1, eq2], and dimensions x = [x1, x2, ...]

function [i, x_now] = findFeasableIntercept(e, c, i, x)

x_now = zeros(size(x));
found = false;

while ~found
    S = solve(e{i(1)}(x), e{i(2)}(x));
    
    % Flexible with dimensionality assuming all variables are named x1, x2, x3, etc...
    if ~isempty(fieldnames(S)) % Account for possibility of lines that don't intercept
        for j = 1:length(x)
            index = "S.x" + sprintf("%d",j);
            x_now(j) = eval(index);
        end
   
        found = isFeasable(x_now, c);
    else
        found = 0;
    end
    
    if found
        break
    end
    
    % Cycle through all combinations of intercepting lines/equations
    i = iCycle(i, e);
    if (class(i) == "string") && (i == "cycle_end")
        sprintf("Feasable region does not exist")
        break
    end
    
end