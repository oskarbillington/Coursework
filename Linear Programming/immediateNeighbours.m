% Find the immediately neighbouring interception points from a list of
% potential ones, by looking along both currently intercepting lines and
% choosing the closest intercept in each direction based on Euclidian
% norms. Input potential neighbours{line1 = [line1, other line, x1, x2,
% ...; repeat], line2 = [-:-]}, current position (i, x_now), and dimensions
% [x1, x2, x3, ...]. Output vector with immediately neighbouring 
% interception points and the distance to the intercept [line a, line b, 
% distance; ...].


function nv_i = immediateNeighbours(nv, i, x_now, x)

min_distance = ones(4,3)*inf; % [line1, other line+, min_distance line1+; line1, other line-, dist; & for line2.]
min_distance([1 2],1) = i(1);
min_distance([3 4],1) = i(2);

for a = 1:2 % For each current line
    
    direction = zeros(size(x));
    
    for k = 1:size(nv{a},1) % For each of the intercepts with other lines
        distance = norm((nv{a}(k,[3,end])) - x_now); % Find closest point in euclidian sense:
        if (distance < min_distance(2*(a-1)+1,3)) % @@@@ Aaah let's assume no 3 lines intercept at one point! @@@@
            min_distance(2*(a-1)+1,3) = distance;
            min_distance(2*(a-1)+1,2) = nv{a}(k,2);
            direction = ((nv{a}(k,[3,end])) - x_now)/distance; % Normal direction of currently found minimum distance
        end
    end
    
    % Now find closest point in the other direction
    
    for k = 1:size(nv{a},1)
        if nv{a}(k,2) ~= min_distance(2*(a-1)+1,2) % Only checking lines not already used @@ unnecessary if already checking directions, except for when the distance was 0! maybe this will help with triple-points...
            distance = norm((nv{a}(k,[3,end])) - x_now);
            new_direction = ((nv{a}(k,[3,end])) - x_now)/distance;
            if (new_direction == direction*(-1)) % Now checking that direction is opposite
                if (distance < min_distance(2*(a-1)+2,3))
                    min_distance(2*(a-1)+2,3) = distance;
                    min_distance(2*(a-1)+2,2) = nv{a}(k,2);
                end
            end
        end
    end
    
end

% Remove Inf-distances (where there are fewer than 4 neighbours)
%flags = zeros(1, size(min_distance,1));

for k = flip(1:size(min_distance,1))
    if (min_distance(k,3) == Inf)
        min_distance(k,:) = [];
    end
end

nv_i = min_distance;

end