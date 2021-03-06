% 1D* LINEAR ELEMENTS**

% (Animate over some variable)
for my_number = 1:15
    
    
    % Environment
    L = 1; % Length*
    Start = 0; % Base coordinate*
    A = 1e-2; % Cross-sectional area
    E = 210e9; % Young's Modulus
    b = 7.85e3*9.81*A; % Distributed axial force (mild steel under self-weight)
    P = 10*9.81; % End applied axial load
    %P = my_number*9.81;
    %N1 = 4; % Number of elements
    N1 = my_number;
    
    n=2e2; % Number of test points for error calculation
    
    Boundaries = [Start, Start+L];
    N0 = N1 + 1; % Number of nodes
    h = L/N1; % Element size
    
    % Dirichlet = {[1, 0], [N0, 0]}; % Dirichlet BCs [node, value]
    Dirichlet = {[1, 0]};
    Neumann = {[N0, P/(A*E)]}; % Neumann BCs [node, value]
    
    gnodes = zeros(N0,1); % List: Global nodes (1 coordinate per node*)
    elements = zeros(N1,2); % List: Element nodes (2 nodes per element**)
    
    for i = 1:size(gnodes,1) % Fill in lists:
        gnodes(i,1) = Boundaries(1) + (i-1)*h;
    end
    for i = 1:size(elements,1)
        elements(i,1) = i;
        elements(i,2) = i+1;
    end
    
    
    % Analytical Solution:
    Solution = @(x) b*L^2/(A*E)*(x/L)*(1-x/(2*L))+P*x/(A*E);
    
    
    % Initiate Stiffness (K) and Force (F) vectors:
    K = zeros(N0);
    F = zeros(N0,1);
    
    
    % Calculate local stiffness matrix and nodal force vectors (same for all
    % elements)
    K_local = zeros(2,2);
    Gradients = {-1/h; 1/h}; % Gradients of local shape functions (trivial: parent elements & integration scheme unnecessary)
    
    for i = 1:size(K_local,1)
        for j = 1:size(K_local,2)
            K_local(i,j) = Gradients{i}*Gradients{j}*A*E*h;
        end
    end
    
    F_local = ones(2,1)*b*h/2;
    
    
    % Assemble
    for e = 1:N1
        ij = elements(e,:);
        
        % Stiffness matrix:
        for local_i = 1:size(K_local,1)
            for local_j = 1:size(K_local,2)
                K(ij(local_i),ij(local_j)) = K(ij(local_i),ij(local_j)) + K_local(local_i, local_j);
            end
        end
        
        % Nodal force vector:
        for k = 1:length(ij)
            F(ij(k)) = F(ij(k)) + F_local(k);
        end
    end
    
    
    % Apply Neumann boundary conditions
    for cond = 1:length(Neumann)
        F(Neumann{cond}(1)) = F(Neumann{cond}(1)) + Neumann{cond}(2)*A*E;
    end
    
    
    % Apply Dirichlet boundary conditions
    for cond = 1:length(Dirichlet)
        F = F - K(:,Dirichlet{cond}(1)).*Dirichlet{cond}(2);
    end
    
    for cond = 1:length(Dirichlet)
        K(:,Dirichlet{cond}(1)) = 0;
        K(Dirichlet{cond}(1),:) = 0;
        K(Dirichlet{cond}(1),Dirichlet{cond}(1)) = 1;
        F(Dirichlet{cond}(1)) = Dirichlet{cond}(2);
    end
    
    
    % Solve set of linear equations
    U = K\F;
    
    
    % Calculate stress and strain in each element
    gnodes_displaced = gnodes + U;
    sigma = zeros(N1,1);
    epsilon = zeros(N1,1);
    
    for e = 1:N1
        epsilon(e) = (gnodes_displaced(e+1)-gnodes_displaced(e))/(gnodes(e+1)-gnodes(e)) - 1;
        sigma(e) = E*epsilon(e);
    end
    
    
    % Plotting
    clf
    hold on
    
    fplot(Solution, Boundaries, "Color", "green", "LineWidth", 1.5)
    plot(gnodes, U, 'o-', "Color", "red")
    legend("Analytical solution", "Approximate function", "Location", "northwest")
    title("Finite element approach")
    %axis([Start, Start+L, -1e-6, 1.2e-6])
    
    y_approx = zeros(1,n); % Average error (%):
    y_solution = zeros(1,n);
    error_vec = zeros(1,n);
    x_now = 0;
    x_step = L/(n-1);
    
    for k = 1:n
        y_approx(k) = Approximate_function(x_now, gnodes, elements, U, h, N1);
        y_solution(k) = Solution(x_now);
        error_vec(k) = y_approx(k) - y_solution(k);
        x_now = x_now + x_step;
    end
    %plot(linspace(Start, Start+L, n), y_approx, 'o', "Color", "black")
    
    error_avg_percent = mean(abs(error_vec./(y_solution+1e-15).*100));
    info1 = sprintf("Mean error: %3.3f %%", error_avg_percent);
    position1 = [min(xlim) max(ylim)]+[diff(xlim)*0.05 -diff(ylim)*0.28];
    text(position1(1), position1(2), info1)
    info2 = sprintf("# of elements: %d", N1);
    position2 = [min(xlim) max(ylim)]+[diff(xlim)*0.05 -diff(ylim)*0.20];
    text(position2(1), position2(2), info2)
    
    drawnow
    pause(0.6)
    
end

% Assemble piece-wise linear Approximate function
function y = Approximate_function(x, gnodes, elements, U, h, N1)

for e = 1:N1
    if ((x-1e-12) <= gnodes(e+1))
        ij = elements(e,:);
        y = U(ij(1)) + (x - h*(e-1))*(U(ij(2)) - U(ij(1)))/h; % Linear approx
        break
    end
end

end
