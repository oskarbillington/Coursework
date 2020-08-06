% Solving a linear programming problem using the SIMPLEX METHOD (by choice)


% Setup
syms x1 x2 real
x = [x1, x2];
max_iter = 4; % Maximum search iterations

% Constraints
c = {@(x) (2*x(1) + 8*x(2) <= 60); @(x) (5*x(1) + 2*x(2) <= 60); @(x) (x(1) >= 0); @(x) (x(2) >= 0); @(x) (x(1) <= 11)}; % c{1} should be <= 60!

% Equations equivalent
e = {@(x) (2*x(1) + 8*x(2) == 60); @(x) (5*x(1) + 2*x(2) == 60); @(x) (x(1) == 0); @(x) (x(2) == 0); @(x) (x(1) == 11)}; % c{1} should be <= 60!

% Task: max f(x1, x2)
f = @(x) 40*x(1) + 88*x(2);



% ALGORITHM: 1) Find a feasable solution by checking if random vertices (found by solving
% random e{i}, e{j}) solve all constraints c{i}. Save the vertex' value to
% f_now and it's coordinates to x_now;

% Initial equation indices (always size 2 independent of size(x))
i = [4,5];

% Cycle through intercepts until one is found in the feasable region
[i, x_now] = findFeasableIntercept(e, c, i, x);

f_now = f(x_now);

% Initialise vector to store iterations: @@@ Assuming optimum function
% value in feasable region > 0
history = zeros(max_iter, 1+length(x)); % [value, x1, x2, ...]

% Initialise variables to store data about next vertex in search path
x_next = zeros(size(x));
i_next = [0 0];
f_next = 0; % @@@ What if zero is greater than the greatest function value in the feasable region?

for iter = 1:max_iter
    
    history(iter,1) = f_now;
    history(iter,[2,end]) = x_now;
    
    % ALGORITHM: 2) Find neighbour vertices (@@@@ HOW!? @@@@) (let's try my_method 2, see paperwork)
    
    % Find all intercepts between lines intercepting at the current point and all other lines
    nv = findIntercepts(e, i, x);
    
    % Find immediately neighbouring vertices @@@ Assuming no triple-points! @@@
    nv_i = immediateNeighbours(nv, i, x_now, x);
    
    
    % ALGORITHM: 3) Initialise f_next = 0 and x_next = empty. (@@@ done outside main loop) Loop for each neighbour vertex: if it is feasable,
    % calculate it's cost f(vertex(i)); if > f_next, f_next = f(vertex(i)) &
    % x_next = x(i). After looping through neighbours, f_next & x_next describe the greatest neighbouring value.
    
    for k = 1:size(nv_i,1)
        S = solve(e{nv_i(k,1)}(x), e{nv_i(k,2)}(x));
        x_check = zeros(size(x));
        
        for j = 1:length(x) % For each dimension
            index = "S.x" + sprintf("%d",j);
            x_check(j) = eval(index); % Find intercept coordinate
        end
        
        if isFeasable(x_check, c) % If the intercept is in the feasable region
            val_check = f(x_check);
            if val_check > f_next
                f_next = val_check;
                i_next = nv_i(k,[1 2]);
                x_next = x_check;
            end
        end
        
    end
    
    
    % ALGORITHM: 4) If f_next > f_now, x_now = x_next and repeat. Else, f_now is the
    % global maximum and the problem is solved.
    
    if f_next > f_now
        f_now = f_next;
        i = i_next;
        x_now = x_next;
    else
        sprintf("Optimum value found: %.2f, at coordinates: %.2f, %.2f", f_now, x_now(1), x_now(2))
        break
    end
    
end


% Plot results @@@ Add better plotting showing values and step directions

while history(end,1) == 0
    history(end,:) = [];
end

my_colours = hot;
plot(history(:,2), history(:,3), "o-", "Color", my_colours(1,:))
axis([-1,20,-1,15])
drawnow

% Or check surfc(history) ...
