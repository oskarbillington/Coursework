% Input phi (yaw), theta (forward/backward wheel roll), and psi (pitch)
% logs from simulation (where column 1 is the time stamps, column 2 the
% corresponding values, and column 3 optionally contains the reference
% values). Outputs animated robot in 3D (wheels w/axle and line from axle
% to *center* of robot body), and the path in 2D

function display = animate_results_3D(log_phi, log_theta, log_psi)

%################################
% User-specified parameters:
Arena_size = 2; % Scales the arena size
Animation_speed = 1; % Ratio of simulation-time to animation-time
Computer_speed = 2.5; % Lower if animation-time runs inconsistently. Increase to plot more details. Controls a filter that can skip close-in-time data points
Camera_follow_robot = false; % Set to 'true' to center the camera on the robot. Set to 'false' to see the whole arena
Grid_on = false; % Toggle 3D arena grid. When turned on, ground colour gradient is turned off (the workaround to be able to stack them would noticably slow down the animation)
Ghost3D = true; % Set to 'true' if the inputs have a 3rd column for the reference signal; the reference-ghost will be plotted with the 3D simulation animation. Otherwise set to 'false'
Ghost2D = false; % ...set to 'true' to also plot the reference path on the 2D plot. Otherwise set to 'false'
%################################

% Setup
time = log_phi(:,1); % == log_theta(:,1) == log_psi(:,1)
phi = log_phi(:,2); % Assuming second column of inputs is simulation performance, not the reference signal (if both were fed into the scope)
theta = log_theta(:,2);
psi = log_psi(:,2);
if Ghost3D == true || Ghost2D == true
    ghost_phi = log_phi(:,3);
    ghost_theta = log_theta(:,3);
    ghost_psi = log_psi(:,3);
end

n = length(time);

my_colour = hot;
delta_colour = length(my_colour)/n; % Because number of colours doesn't match n

R = 0.027; % Wheel radius
W = 0.14; % Body width
H = 0.22; % Body height

x_m = 0; % Start at origin
y_m = 0;
z_m = R;

if Ghost3D == true || Ghost2D == true
    ghost_x_m = 0; % Start at origin
    ghost_y_m = 0;
    ghost_z_m = R;
end

d_max = Arena_size*R*max(theta); % Display boundaries

% Initiate figure, time controls, and draw arena-boundary
display = figure;
time1 = 0;
tic
first_plot = 1;

subplot(2,1,2);
hold on
patch([d_max -d_max -d_max d_max], [d_max d_max -d_max -d_max], [0.131, 0.139, 0.139])
patch(0.99*[d_max -d_max -d_max d_max], 0.99*[d_max d_max -d_max -d_max], [0.9, 0.9, 0.9])

for i = 2:n
    
    if time1 == 0 % Time of last plotted frame
        time1 = time(i-1);
    end
    
    if ((time(i)-time1) > 1e-1*Animation_speed/Computer_speed) || i == n || i == 2 % If enough simulation time has passed since last plotted frame
        
        % Plot 3D model
        h1 = subplot(2,1,1);
        cla(h1)
        hold on
        
        if Ghost3D == true && first_plot ~= 1 % Plot ghost underneath:
            
            ghost_w_x = [ghost_x_m + W/2*sin(ghost_phi(i)), ghost_x_m - W/2*sin(ghost_phi(i))]; % Wheel axis-line positions:
            ghost_w_y = [ghost_y_m - W/2*cos(ghost_phi(i)), ghost_y_m + W/2*cos(ghost_phi(i))];
            ghost_w_z = [ghost_z_m, ghost_z_m];
            ghost_n_wheel = [ghost_w_x(1)-ghost_w_x(2), ghost_w_y(1)-ghost_w_y(2), 0]; % Wheel circles normal vector
            ghost_b_x = [ghost_x_m, ghost_x_m + H/2*sin(ghost_psi(i))*cos(ghost_phi(i))]; % Body center-line positions:
            ghost_b_y = [ghost_y_m, ghost_y_m + H/2*sin(ghost_psi(i))*sin(ghost_phi(i))];
            ghost_b_z = [ghost_z_m, ghost_z_m + H/2*cos(ghost_psi(i))];
            ghost_center_r = [ghost_w_x(1), ghost_w_y(1), ghost_w_z(1)];
            ghost_center_l = [ghost_w_x(2), ghost_w_y(2), ghost_w_z(2)];
            
            plot3(ghost_w_x,ghost_w_y,ghost_w_z, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 5) % Draw robot:
            plot3(ghost_b_x,ghost_b_y,ghost_b_z, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 5)
            
            ghost_angle=0:0.01:2*pi; % Draw wheels:
            ghost_v=null(ghost_n_wheel);
            ghost_points_r=repmat(ghost_center_r',1,size(ghost_angle,2))+R*(ghost_v(:,1)*cos(ghost_angle)+ghost_v(:,2)*sin(ghost_angle));
            ghost_points_l=repmat(ghost_center_l',1,size(ghost_angle,2))+R*(ghost_v(:,1)*cos(ghost_angle)+ghost_v(:,2)*sin(ghost_angle));
            plot3(ghost_points_l(1,:),ghost_points_l(2,:),ghost_points_l(3,:),'-','Color',[0.9, 0.9, 0.9],'LineWidth',3);
            plot3(ghost_points_r(1,:),ghost_points_r(2,:),ghost_points_r(3,:),'-','Color',[0.9, 0.9, 0.9],'LineWidth',3);
            
            ghost_l_x = ghost_w_x + R*cos(ghost_phi(i))*sin(ghost_theta(i)); % Draw line on wheels
            ghost_l_y = ghost_w_y + R*sin(ghost_phi(i))*sin(ghost_theta(i));
            ghost_l_z = ghost_w_z + R*cos(ghost_theta(i));
            plot3([ghost_w_x(1), ghost_l_x(1)],[ghost_w_y(1), ghost_l_y(1)],[ghost_w_z(1), ghost_l_z(1)], '-','Color',[0.9, 0.9, 0.9],'LineWidth',3);
            plot3([ghost_w_x(2), ghost_l_x(2)],[ghost_w_y(2), ghost_l_y(2)],[ghost_w_z(2), ghost_l_z(2)], '-','Color',[0.9, 0.9, 0.9],'LineWidth',3);
            
        end
        
        w_x = [x_m + W/2*sin(phi(i)), x_m - W/2*sin(phi(i))]; % Wheel axis-line positions:
        w_y = [y_m - W/2*cos(phi(i)), y_m + W/2*cos(phi(i))];
        w_z = [z_m, z_m];
        n_wheel = [w_x(1)-w_x(2), w_y(1)-w_y(2), 0]; % Wheel circles normal vector
        b_x = [x_m, x_m + H/2*sin(psi(i))*cos(phi(i))]; % Body center-line positions:
        b_y = [y_m, y_m + H/2*sin(psi(i))*sin(phi(i))];
        b_z = [z_m, z_m + H/2*cos(psi(i))];
        center_r = [w_x(1), w_y(1), w_z(1)];
        center_l = [w_x(2), w_y(2), w_z(2)];
        
        plot3(w_x,w_y,w_z, 'Color', 'black', 'LineWidth', 5) % Draw robot:
        plot3(b_x,b_y,b_z, 'Color', 'blue', 'LineWidth', 5)
        
        angle=0:0.01:2*pi; % Draw wheels:
        v=null(n_wheel);
        points_r=repmat(center_r',1,size(angle,2))+R*(v(:,1)*cos(angle)+v(:,2)*sin(angle));
        points_l=repmat(center_l',1,size(angle,2))+R*(v(:,1)*cos(angle)+v(:,2)*sin(angle));
        plot3(points_l(1,:),points_l(2,:),points_l(3,:),'r-','LineWidth',3);
        plot3(points_r(1,:),points_r(2,:),points_r(3,:),'r-','LineWidth',3);
        
        l_x = w_x + R*cos(phi(i))*sin(theta(i)); % Draw line on wheels
        l_y = w_y + R*sin(phi(i))*sin(theta(i));
        l_z = w_z + R*cos(theta(i));
        plot3([w_x(1), l_x(1)],[w_y(1), l_y(1)],[w_z(1), l_z(1)], 'r-','LineWidth',3);
        plot3([w_x(2), l_x(2)],[w_y(2), l_y(2)],[w_z(2), l_z(2)], 'r-','LineWidth',3);
        
        title(sprintf('ROBOT SIMULATION\nTime: %.2f/%d sec,  Pitch: %.2f deg', time(i), floor(time(end)), psi(i)/pi*180))
        
        if Grid_on == true
            grid on
        else
            patch([d_max -d_max -d_max d_max], [d_max d_max -d_max -d_max], [0 0 0 0], [1 1 -1 -1])
        end
        
        if Camera_follow_robot == true
            axis([x_m-3*W, x_m+3*W, y_m-3*W, y_m+3*W, 0, H*0.6])
        else
            axis([-d_max d_max -d_max d_max 0 H])
        end
        view(45,30)
        drawnow limitrate
        
        % Plot 2D path
        subplot(2,1,2);
        hold on
        if Ghost2D == true
            plot(ghost_x_m,ghost_y_m,'.', 'Color', [1,1,1],'LineWidth', 2)
        end
        plot(x_m,y_m,'.', 'Color', my_colour(ceil(i*delta_colour),:),'LineWidth', 2)
        axis([-3.5*d_max 3.5*d_max -1.05*d_max 1.05*d_max])
        xlabel(sprintf('x_m\nUnits: meters'))
        ylabel('y_m')
        title('Path history','LineWidth',0.3,'FontSize',10)
        drawnow limitrate
        
        % Animation speed to scale (assuming really good computer speed)
        t = toc;
        pause((time(i)-time1-t)/Animation_speed)
        
        % Reset time-controls
        tic
        time1 = 0;
        
        % On first iteration, ask user to start animation
        if first_plot == 1
            subplot(2,1,1)
            if Camera_follow_robot == true
                textbox = text(0,0,H*0.5,'Press a key to start animation','Color','green','FontSize',15);
            else
                textbox = text(0,0,H*0.9,'Press a key to start animation','Color','green','FontSize',15);
            end
            set(textbox,'visible','on','HorizontalAlignment','center','VerticalAlignment','middle')
            pause
            delete(textbox)
            first_plot = 0;
        end
        
    end
    
    % Update x_m, y_m
    x_m = x_m + R*(theta(i)-theta(i-1))*cos(phi(i));
    y_m = y_m + R*(theta(i)-theta(i-1))*sin(phi(i));
    
    if Ghost3D == true || Ghost2D == true
        ghost_x_m = ghost_x_m + R*(ghost_theta(i)-ghost_theta(i-1))*cos(ghost_phi(i));
        ghost_y_m = ghost_y_m + R*(ghost_theta(i)-ghost_theta(i-1))*sin(ghost_phi(i));
    end
    
    % Failure if the robot falls over
    if abs(psi(i)) >= pi/2
        break
    end
    
end

end