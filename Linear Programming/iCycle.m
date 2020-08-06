% Cycles through all combinations of (i, j~=i) in ascending order, and finally returns
% "cycle_end". Input current i = [int1, int2] and e = {equations to intercept}

function i = iCycle(i, e)

if i(2) < length(e)
    i(2) = i(2) + 1;
elseif i(1) < (length(e)-1)
    i(1) = i(1) + 1;
    i(2) = i(1) + 1;
else
    i = "cycle_end";
end

end