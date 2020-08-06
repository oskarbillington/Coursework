% 1D* QUADRATIC ELEMENTS**


% (Animate over some variable)
for my_number = 0:10
    
    
    % Environment
    L = 1; % Length*
    Start = 0; % Base coordinate*
    A = 1e-2; % Cross-sectional area
    E = 210e9; % Young's Modulus
    b = 7.85e3*9.81*A; % Distributed axial force (mild steel under self-weight)
    %P = 10*9.81; % End applied axial load
    P = my_number*9.81;
    N1 = 5; % Number of elements
    %N1 = my_number;
    
    n = 1e5; % Number of error-calculation points
    
    Boundaries = [Start, Start+L];
    N0 = 2*N1 + 1; % Number of nodes
    h = L/N1; % Element size
    d = h/2; % Distance between nodes
    
    %Dirichlet = {[1, 0], [N0, 0]}; % Dirichlet BCs [node, value]
    %Neumann = {};
    Dirichlet = {[1, 0]};
    Neumann = {[N0, P/(A*E)]}; % Neumann BCs [node, value]
    
    gnodes = zeros(N0,1); % List: Global nodes (1 coordinate per node*)
    elements = zeros(N1,3); % List: Element nodes (4 nodes per element**) 
    
    for i = 1:size(gnodes,1) % Fill in lists:
        gnodes(i,1) = Boundaries(1) + (i-1)*d;
    end
    for i = 1:size(elements,1)
        elements(i,1) = 1 + 2*(i-1);
        elements(i,2) = elements(i,1) + 2;
        elements(i,3) = elements(i,1) + 1;
    end
    
    
    
    % Analytical Solution:
    Solution = @(x) b*L^2/(A*E)*(x/L)*(1-x/(2*L))+P*x/(A*E);
    
    
    % Calculate local stiffness matrix and nodal force vectors (same for all
    % elements)
    K_local = zeros(3,3);
    F_local = zeros(3,1);
     
    syms x z real
    Local_shape_funcs = {2*x^2/h^2-3*x/h+1, 2*x^2/h^2-x/h, -4*x^2/h^2+4*x/h}; % (Local coordinates: [1-3-2])
    Local_gradients = Local_shape_funcs; % (Initialise) local gradients
    for i = 1:length(Local_gradients)
        Local_gradients{i} = diff(Local_shape_funcs{i});
        Local_shape_funcs{i} = subs(Local_shape_funcs{i}, x, (1+z)*h/2);
        Local_gradients{i} = subs(Local_gradients{i}, x, (1+z)*h/2);
    end
    
    Je = h/2; % Coordinate transformation; Jacobian determinant
    
    for i = 1:size(K_local,1)
        for j = 1:size(K_local,2)
            K_local(i,j) = A*E*Gauss_Legendre(Local_gradients{i}*Local_gradients{j}*Je);
        end
    end
    
    for i = 1:length(F_local)
        F_local(i) = b*Gauss_Legendre(Local_shape_funcs{i}*Je);
    end
    
    
    % Assemble
    K = zeros(N0);
    F = zeros(N0,1);
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
    if ~isempty(Neumann)
        for cond = 1:length(Neumann)
            F(Neumann{cond}(1)) = F(Neumann{cond}(1)) + Neumann{cond}(2)*A*E;
        end
    end
    
    
    % Apply Dirichlet boundary conditions
    if ~isempty(Dirichlet)
        for cond = 1:length(Dirichlet)
            F = F - K(:,Dirichlet{cond}(1)).*Dirichlet{cond}(2);
        end
        for cond = 1:length(Dirichlet)
            K(:,Dirichlet{cond}(1)) = 0;
            K(Dirichlet{cond}(1),:) = 0;
            K(Dirichlet{cond}(1),Dirichlet{cond}(1)) = 1;
            F(Dirichlet{cond}(1)) = Dirichlet{cond}(2);
        end
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
    fplot(Solution, Boundaries, "Color", "yellow", "LineWidth", 2.3)
    plot(gnodes, U, 'o', "Color", "red") %@@@ need way to plot nicely
    title("Quadratic 1D finite elements")
    %axis([Start, Start+L, -1e-6, 1.2e-6])
    
    my_colour = cool; % Plotting the element piece-wise quadratic polynomial approximations:
    gr = floor(size(my_colour,1)/N1);
    y_approx = cell(N1);
    for e = 1:N1
        ijk = [1, 2, 3] + 2*(e-1);
        P = polyfit(gnodes(ijk), U(ijk), 2);
        y_approx{e} = @(x) P(1)*x^2 + P(2)*x + P(3);
        fplot(y_approx{e}, [gnodes(ijk(1)), gnodes(ijk(3))], "Color", my_colour(gr*e,:))
    end
    
    legend("Analytical solution", "Approximate nodes", "Location", "northwest")
    
    info2 = sprintf("# of elements: %d", N1);
    position2 = [min(xlim) max(ylim)]+[diff(xlim)*0.05 -diff(ylim)*0.20];
    text(position2(1), position2(2), info2)
    
    drawnow
    pause(0.3)
    
end

% Integrates upto 3rd order polynomial exactly from -1 to 1. Input integrand output integral. Gauss-Legendre method.
function integral = Gauss_Legendre(integrand)

integral = 0;
k = [0.577350269189626];
w = [1];

for i = 1:length(k)
    z = k;
    integral = integral + w(i)*subs(integrand);
    z = -k;
    integral = integral + w(i)*subs(integrand);
end

end

