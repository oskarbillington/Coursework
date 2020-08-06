% Input list of equations e{}, current location as indeces of currently
% intercepting lines i=[line1, line2], and dimension variables x = [x1, x2,
% ...]. Outputs lists of interception points for both lines with every
% other line [some current line, some other line, x1, x2, x3, ...]. These
% points represent potential nearest-neighbours.

function nv = findIntercepts(e, i, x)

nv = {zeros(length(e)-2, length(x)+2), zeros(length(e)-2, length(x)+2)}; % Potential neighbouring vertices: [line1, line2, x1, x2, x3, ...]
k = 1;

for j = 1:length(e)
    if ~ismember(i,j) % Looking at all lines than the currently intersecting ones
        
        for a = 1:2 % For each line at the current intercept
            nv{a}(k,[1,2]) = [i(a),j]; % Save new line indeces
            S = solve(e{i(a)}(x), e{j}(x)); % ... and new corresponding intercept coordinates
            if ~isempty(fieldnames(S)) % Account for possibility of lines that don't intercept
                for b = 1:length(x) % Then for each dimension:
                    index = "S.x" + sprintf("%d",b);
                    nv{a}(k,2+b) = eval(index); % Save coordinate
                end
            else
                nv{a}(k,:) = [];
            end
        end
        k = k+1;
        
    end
end

end