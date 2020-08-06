% LEGO Football CWM UOXF -- Trinity Term 2020


% @@@ Initialisation @@@
g = 9.81;						% gravity acceleration [m/sec^2]

% Physical Parameters
m = 0.024;  					% wheel weight [kg]
R = 0.027;  					% wheel radius [m]
Jw = m*R^2/2;                   % wheel inertia moment [kgm^2]
M = 0.6;    					% body weight [kg]
W = 0.14;   					% body width [m]
D = 0.04;   					% body depth [m] 0.04 changed
H = 0.22;                   	% body height [m]
L = H/2;    					% distance of the center of mass from the wheel axle [m]
Jpsi = M*L^2/3;                 % body pitch inertia moment [kgm^2]
Jphi = M*(W^2+D^2)/12;          % body yaw inertia moment [kgm^2]
b = 0.0022;                    	% friction coefficient between body & DC motor [Nms/rad]
b_bar = 0.0022; 			    % friction coefficient due to rotation along the wheel axis [Nms/rad]

% DC Motor Parameters (assumed to be similar to NXT after comparison experiments)
Rdc = 6.69; 					% DC motor resistance [ƒ¶]
Ke = 0.468;                   	% DC motor back EMF constant [Vsec/rad]
Kt = 0.317;                   	% DC motor torque constant [Nm/A]

% Other
plant_dim_in = 2; % Number of plant inputs
plant_dim_out = 3;


% @@@ Linear System Matrices (derived by hand from physical system modelling) @@@
q = ((2*m + M)*R^2 + 2*Jw)*(M*L^2 + Jpsi) - M^2*R^2*L^2;

A1 = [0, 0, 1, 0;...
      0, 0, 0, 1;...
      0, -g*R*M^2*L^2/q, -2*(M*R*L + M*L^2 + Jpsi)*(Kt*Ke/Rdc + b)/q - 2*(M*L^2 + Jpsi)*b_bar/q, 2*(M*R*L + M*L^2 + Jpsi)*(Kt*Ke/Rdc + b)/q;...
      0, g*M*L*((2*m + M)*R^2 + 2*Jw)/q, 2*((2*m + M)*R^2 + 2*Jw + M*R*L)*(Kt*Ke/Rdc + b)/q + 2*M*R*L*b_bar/q, -2*((2*m + M)*R^2 + 2*Jw + M*R*L)*(Kt*Ke/Rdc + b)/q];

B1 = [0; 0; 2*(M*R*L + M*L^2 + Jpsi)*Kt/Rdc/q; -2*((2*m + M)*R^2 + 2*Jw + M*R*L)*Kt/Rdc/q];

C1 = [1, 0, 0, 0;...
      0, 0, 0, 1];
  
D1 = [0;0];

plant1 = ss(A1,B1,C1,D1);

A2 = [0, 1;...
      0, -W^2/2/R^2*(Kt*Ke/Rdc + b + b_bar)/(1/2*m*W^2 + Jphi + Jw*W^2/2/R^2)];

B2 = [0; -W/R*Kt/Rdc/(1/2*m*W^2 + Jphi + Jw*W^2/2/R^2)];

C2 = [1 0];

D2 = 0;

plant2 = ss(A2,B2,C2,D2);

% (Unnecessary:) Assembly. State x = [theta, psi, dot_theta, dot_psi, phi, dot_phi]', input u = [u1, u2]', ouput y = [theta, dot_psi, phi]'.
%A = [A1, zeros(size(A1,1), size(A2,2));...
%     zeros(size(A2,1), size(A1,2)), A2];
%B = [B1, zeros(size(B1,1), size(B2,2));...
%     zeros(size(B2,1), size(B1,2)), B2];
%C = [C1, zeros(size(C1,1), size(C2,2));...
%     zeros(size(C2,1), size(C1,2)), C2];
%D = zeros(plant_dim_out, plant_dim_in); 
%Plant = ss(A,B,C,D);
  

% @@@ Control design @@@

% dot_psi PI-controller
z_PI = 7; % = K_I_dot_psi/K_P_dot_psi. Set to be smaller than the least negative pole of G_dot_psi
tfs = zpk(plant1);
G_dot_psi = tfs(2);
s = tf('s');
%L = minreal(-G_dot_psi*(s+z_PI)/s); % Open loop transfer function with unity K_P-gain
%figure
%rlocus(L)
K_P_dot_psi = -50; % < -11.8 to ensure (marginal) stability. Larger value gives better GM; PM converges from 106deg at marginal stability to 90deg as abs(K)->inf
C_dot_psi = K_P_dot_psi*(s + z_PI)/s; % PI-controller
%L_dot_psi = minreal(C_dot_psi*G_dot_psi); % Open-loop PI transfer function
%[GM_PI,PM_PI] = margin(L_dot_psi)

% New virtual plant including the PI-controller:
G = [tfs(1); tfs(2)];
G_tilde = minreal([1 0] * G / (1 + C_dot_psi*[0 1]*G));

% theta PID-controller
w0 = 1; % PID zeros found by projecting performance specs onto natural frequency and damping ratio:
xi = 0.8;
r_ID = w0^2;
r_PD = 2*xi*w0; 
%L = minreal(G_tilde*(s^2 + r_PD*s + r_ID)/s);
%figure
%rlocus(-L)
K_D_theta = -8; % > -10.5 to ensure stability.
C_theta = K_D_theta*(s^2 + r_PD*s + r_ID)/s;
%[GM_PID,PM_PID] = margin(minreal(C_theta*G_tilde))

% Full transfer function theta/theta_ref:
TF_theta = minreal(C_theta*G_tilde/(1 + C_theta*G_tilde)); % Negative feedback
%figure
%step(TF_theta)
%hold on
%fplot(@(x) 1.05, 'Color', 'green')
%fplot(@(x) 0.95, 'Color', 'green')
%figure
%nyquist(minreal(C_theta*G_tilde))
%grid on

% Yaw proportional controller
G_phi = zpk(plant2);
%figure
%rlocus(minreal(-G_phi))
K_phi_OL = -A2(2,2)/B2(2); % Open-loop gain
K_P_phi = -50; % This found by root-locus (found any negative number works)
L = G_phi*K_P_phi;
%figure
%margin(L);
TF_phi_OL = K_phi_OL*G_phi;
TF_phi = minreal(K_P_phi*G_phi/(1+K_P_phi*G_phi));

