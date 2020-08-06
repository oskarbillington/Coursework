% 1D* LINEAR ELEMENTS**


% Environment 
L = 1; % Length*
Start = 0; % Base coordinate*
Boundaries = [Start, Start+L];
b = 1; % Distributed axial force
A = 1; % Cross-sectional area
E = 1; % Young's Modulus
N1 = 2; % Number of elements
N0 = N1 + 1; % Number of nodes
h = L/N1; % Element size

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
Solution = @(x) L^2*b/(2*A*E)*(x/L)*(1-x/L);  


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


% Dirichlet boundary conditions
Dirichlet = {[1, 0], [N0, 0]}; % [node, value]

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














