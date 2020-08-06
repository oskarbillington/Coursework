% Check if input1 = x_now satisfies each input2 = constraints{} and returns
% boolean

function feasable = isFeasable(x_now, constraints)

feasable = 1;

for i = 1:length(constraints)
    feasable = feasable*isAlways(constraints{i}(x_now));
end
    
end