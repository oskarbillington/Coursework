function SpheroidKicker2000()
% Free kick game. Use left mouse button in the menu and any key to play.

% - Clean up. Combine variables? Comment.

% Known bugs: (1) quitting while shot is in the air (story mode) and
% reentering the game puts you back at that level, allowing "cheating";
% (2) after selecting game mode, the screen must be clicked to be focused;
% (3) entering anything but a positive real number in practise mode yields
% error

global Level difficulty FPS quitting Practise_input Practise_mode distance_global score d key_press menu_press
FPS = 50;
clf

%%% Variables control: execute the game mode chosen in the last iteration %%%
quitting = 0;
key_press = 0;
menu_press = 0;
if score < 999
    Practise_mode = 0;
end
if Practise_mode == 1
    Level = 8;
elseif isempty(score) || score == 0
    Level = 3;
elseif difficulty == 1
    Level = 3+0.5*score;
elseif difficulty == 2
    Level = 3+score;
end

%%% Run Main Menu unless the last shot scored a point %%%
if isempty(score) == 1 || score < 1
    score = 0;
    d = 20; % d = distance from kickoff to target
    
    % Graphics: design vectors
    goal = [d,(d*1.003+0.1),(d*1.003+0.1),d,d;2+0.02*d,2+0.02*d,4+0.04*d,4+0.04*d,2+0.02*d];
    grass = [-0.1*d,d*1.1,d*1.1,-0.1*d,-0.1*d;(d*0.003+0.5),(d*0.003+0.5),-1-0.1*d,-1-0.1*d,(d*0.003+0.5)];
    xrugby = [4,5,6,6,7,7,8,8,7,7,6,6,5,4,3,3,2,2,1,1,2,2,3,3,4]-4.5;
    yrugby = [1,1,2,3,4,5,6,10,11,12,13,14,15,15,14,13,12,11,10,6,5,4,3,2,1];
    rugby = (0.02+d/800).*[xrugby;yrugby];
    
    % Graphics: initial plottings
    hold on
    fill(goal(1,:),goal(2,:), 'b');
    fill(grass(1,:),grass(2,:), 'g');
    fill(rugby(1,:),rugby(2,:),'r');
    axis([-0.5-0.005*d,d*1.005+0.5,-0.5-0.005*d,d*1.005+0.5])
    drawnow
    
    % ui buttons for selecting game mode
    easy_button = uicontrol('Style','pushbutton', ...
        'String', 'Easy', ...
        'Position', [240 220 80 24],...
        'Callback', @Easy_Callback);
    hard_button = uicontrol('Style','pushbutton', ...
        'String', 'Hard', ...
        'Position', [240 190 80 24],...
        'Callback', @Hard_Callback);
    practise_button = uicontrol('Style','pushbutton', ...
        'String', 'Practise', ...
        'Position', [240 160 80 24],...
        'Callback', @Practise_Callback);
    quit_button = uicontrol('Style','pushbutton', ...
        'String', 'Quit', ...
        'Position', [240 130 80 24],...
        'Callback', @Quit_Callback);
    
    % Title animation
    number1 = 1;
    c_max = 24;
    title_colour = cool(c_max);
    TitleBox = annotation('textbox',[.257 .4 .2 .3],...
        'FitBoxToText','on','EdgeColor','w', 'FontWeight', 'Bold',...
        'String','  Spheroid Kicker 2000',...
        'FontSize',18,'Color',title_colour(1,:));
    
    while true
        for c = c_max/4:c_max*3/4
            if menu_press == 1
                break
            end
            set(TitleBox,'Visible','off')
            if number1 == 1
                TitleBox = annotation('textbox',[.257 .4 .2 .3],...
                    'FitBoxToText','on','EdgeColor','w', 'FontWeight', 'Bold',...
                    'String','  Spheroid Kicker 2000',...
                    'FontSize',18,'Color',title_colour(c,:));
            else
                TitleBox = annotation('textbox',[.257 .4 .2 .3],...
                    'FitBoxToText','on','EdgeColor','w', 'FontWeight', 'Bold',...
                    'String','  Spheroid Kicker 2000',...
                    'FontSize',18,'Color',title_colour(c_max-c,:));
            end
            set(TitleBox,'Visible','on')
            pause(0.07)
        end
        number1 = number1*(-1);
        if menu_press == 1
            break
        end
    end
    
    % Resetting the distance d
    if Practise_mode == 0 && difficulty == 1
        d = 5;
        distance_global = d;
    elseif Practise_mode == 0
        d = 10;
        distance_global = d;
    end
    
    % Settings if menu is bypassed
elseif score < 999
    d = distance_global;
else
    score = 0;
end

%%% Responses to button clicks %%%
    function Easy_Callback(~,~)
        difficulty = 1;
        menu_press = 1;
        pause(0.2)
        set(easy_button,'Visible','off')
        set(hard_button,'Visible','off')
        set(practise_button,'Visible','off')
        if Practise_mode == 1
            d = str2double(Practise_input.String);
            set(Practise_input,'Visible','off')
        else
            Practise_mode = 0;
        end
        set(quit_button,'Visible','off')
        set(TitleBox,'Visible','off')
    end

    function Hard_Callback(~,~)
        difficulty = 2;
        menu_press = 1;
        Practise_mode = 0;
        pause(0.2)
        set(easy_button,'Visible','off')
        set(hard_button,'Visible','off')
        set(practise_button,'Visible','off')
        set(quit_button,'Visible','off')
        set(TitleBox,'Visible','off')
    end

    function Practise_Callback(~,~)
        pause(0.2)
        Practise_mode = 1;
        Level = 8;
        set(practise_button,'Visible','off')
        set(hard_button,'Enable','off')
        set(easy_button,'Enable','off')
        Practise_input = uicontrol('Style', 'Edit',...
            'Position', [240 160 80 24],...
            'String', 'Enter distance',...
            'Callback', @Easy_Callback);
    end

    function Quit_Callback(~,~)
        quitting = 1;
        score = 0;
        close
    end

    function Menu_Callback(~,~)
        menu_press = 0;
        score = 0;
        SpheroidKicker2000
    end

    function Retry_Callback(~,~)
        score = 999;
        SpheroidKicker2000
    end

%%% Graphics: game interface %%%
clf
hold on
goal = [d,(d*1.003+0.1),(d*1.003+0.1),d,d;2+0.02*d,2+0.02*d,4+0.04*d,4+0.04*d,2+0.02*d];
grass = [-1-0.1*d,d*1.1+1,d*1.1+1,-0.1*d-1,-0.1*d-1;(d*0.003+0.5),(d*0.003+0.5),-1-0.1*d,-1-0.1*d,(d*0.003+0.5)];
xrugby = [4,5,6,6,7,7,8,8,7,7,6,6,5,4,3,3,2,2,1,1,2,2,3,3,4]-4.5;
yrugby = [1,1,2,3,4,5,6,10,11,12,13,14,15,15,14,13,12,11,10,6,5,4,3,2,1];
rugby = (0.02+d/800).*[xrugby;yrugby];
fill(goal(1,:),goal(2,:), 'b');
fill(grass(1,:),grass(2,:), 'g');
plot2 = fill(rugby(1,:),rugby(2,:),'r');
xarrow = [1,10,10,20,10,10,1,1];
yarrow = [6,6,3,8,13,10,10,6]-8;
arrow = (d^(0.9)/100).*[xarrow;yarrow];

if Practise_mode == 0
    Score_display = annotation('textbox',[.147 .6 .2 .3],...
        'String',sprintf('  Score:  %d    ',score),...
        'FitBoxToText','on','Color','k');
    menu_button = uicontrol('Style','pushbutton', ...
        'String', 'Main Menu', ...
        'Position', [85 325 80 24],...
        'Callback', @Menu_Callback);
else
    menu_button = uicontrol('Style','pushbutton', ...
        'String', 'Main Menu', ...
        'Position', [85 360 80 24],...
        'Callback', @Menu_Callback);
end

%%% Variables and control, arrow animation: angle selection %%%
key_press = 0;
gcf;
set(gcf, 'KeyPressFcn', @key_press_Fcn);
direction = 1;
anglerange = linspace(0.05,0.495*pi,ceil(0.5*FPS/(Level)^(0.5)));
frames_rotate = length(anglerange);
colour1 = autumn(frames_rotate);
MyAngle = [];

%%% Arrow animation: press a key to select the angle %%%
if quitting == 0
    plot3 = fill(arrow(1,:),arrow(2,:),'y');
else
    close
end

while true
    for k1 = 1:(frames_rotate-1)
        if key_press == 1
            if direction == 1
                MyAngle = anglerange(k1);
            elseif direction == -1
                MyAngle = anglerange(frames_rotate-k1);
            end
            pause(0.1)
            break
        end
        set(plot3,'Visible','off')
        arrow = rotateabout(arrow,direction*pi/(2*frames_rotate),0,0);
        colour_k1 = colour1(ceil(frames_rotate*0.5+abs(frames_rotate*0.5-k1)),:);
        plot3 = fill(arrow(1,:),arrow(2,:),colour_k1);
        axis([-0.5-0.005*d,d*1.005+0.5,-0.5-0.005*d,d*1.005+0.5])
        drawnow
        pause(5/((difficulty^5+log(Level))*FPS))
    end
    direction = direction*(-1);
    if key_press == 1 && isempty(MyAngle) == 0
        break
    end
    pause(0.03)
end

%%% Variables and control, arrow animation: force selection %%%
key_press = 0;
direction = 1;
if d <= 40
    powerrange = linspace(-0.2*(difficulty-2)*sqrt(40*9.81*1.5),sqrt(40*9.81*1.5),ceil(FPS*0.65/(1+log(Level))));
else
    powerrange = linspace(-0.2*(difficulty-2)*sqrt(d*9.81*1.5),sqrt(d*9.81*1.5),ceil(FPS*0.65/(1+log(Level))));
end
frames_power = length(powerrange);
colour2 = autumn(frames_power);
MyForce = [];

%%% Arrow animation: initial arrow shrinking %%%
for kx = 1:(0.5*frames_power)
    set(plot3,'Visible','off')
    arrow_scaled = (1-kx/(0.5*frames_power)).*arrow;
    plot3 = fill(arrow_scaled(1,:),arrow_scaled(2,:),colour_k1);
    axis([-0.5-0.005*d,d*1.005+0.5,-0.5-0.005*d,d*1.005+0.5])
    drawnow
    pause(5/((difficulty^5+log(Level))*FPS))
end

%%% Arrow animation: press a key to select the force %%%
while true
    for k2 = 1:frames_power
        if key_press == 1
            if direction == 1
                MyForce = powerrange(k2);
            elseif direction == -1
                MyForce = powerrange(frames_power-k2);
            end
            break
        end
        set(plot3,'Visible','off')
        arrow_scaled = (0.5-0.5*direction+direction*k2/frames_power).*arrow;
        colour_k2 = colour2((frames_power-k2+1)*0.5*(direction+1)+(k2)*(-0.5)*(direction-1),:);
        plot3 = fill(arrow_scaled(1,:),arrow_scaled(2,:),colour_k2);
        axis([-0.5-0.005*d,d*1.005+0.5,-0.5-0.005*d,d*1.005+0.5])
        drawnow
        pause(5/((difficulty^5+log(Level))*FPS))
    end
    direction = direction*(-1);
    if key_press == 1 && isempty(MyForce) == 0
        break
    end
end

%%% Variables and control, calculation and animation: ball path %%%
t = (MyForce)*sin(MyAngle)*(2/9.81); % The time the ball spends in the air
global x d_landing frames_path
d_landing = 2*((MyForce)^(2))*tan(MyAngle)*((cos(MyAngle))^(2))/9.81; % Where the ball lands, function of force and angle
if d_landing > d
    t_goalpass = d/(MyForce*cos(MyAngle)); % The time the ball spends in the displayed area
end
if t >= 2.5 && d < 20 && d_landing >= d % (Attempting a smooth and consistent animation by manipulating the number of frames)
    frames_path = ceil(1*d+5);
    pause_time = t/frames_path;
elseif t < 2.5 && d < 20 && d_landing >= d
    frames_path = ceil(2.5*d);
    pause_time = t/frames_path;
elseif t < 2.5
    if ceil(t/0.05) > 10
        frames_path = ceil(t/0.05);
    else
        frames_path = 10;
    end
    pause_time = 0.05;
else
    frames_path = 50;
    pause_time = t/frames_path;
end

if d_landing > d
    times = linspace(0,t_goalpass,frames_path);
else
    times = linspace(0,t,frames_path);
end
x = times.*(MyForce*cos(MyAngle)); % x-axis positions at realistic times

path = ballpath(MyForce,MyAngle);

%%% x-y-path (neglecting air resistance) from inputs: force (F), angle (a) %%%
    function pathequation = ballpath(F,a)
        g = 9.81;
        if d_landing < d
            pos = linspace(0,1.01*d_landing,frames_path);
        else
            pos = linspace(0,d*1.01,frames_path);
        end
        pathequation = pos.*tan(a)-(pos.^2).*(g/(2*F^2*cos(a)*cos(a)));
    end

%%% Ball animation %%%
no_kick = 0;
if MyForce ~= powerrange(1)
    for k3 = 2:length(x)
        
        if path(k3) <= 0
            set(plot2,'Visible','off')
            rugby = translate(rugby,path(k3-1)*(5*cos(MyAngle)+5)*((x(k3)-x(k3-1))/(path(k3-1)+abs(k3))),-path(k3-1));
            rugby = rotateabout(rugby,nthroot(0.35*pi/2-MyAngle,3)*MyForce*k3/(FPS*(x(k3))^(1.4)),(max(rugby(1,:))+min(rugby(1,:)))/2,(max(rugby(2,:))+min(rugby(2,:)))/2);
            plot2 = fill(rugby(1,:),rugby(2,:),'r');
            axis([-0.5-0.005*d,d*1.005+0.5,-0.5-0.005*d,d*1.005+0.5])
            drawnow
            break
        end
        
        set(plot2,'Visible','off')
        rugby = translate(rugby,(x(k3)-x(k3-1)),(path(k3)-path(k3-1)));
        rugby = rotateabout(rugby,nthroot(0.35*pi/2-MyAngle,3)*MyForce*k3*(0.05*d)/(FPS*(x(k3))^(1.4)),(max(rugby(1,:))+min(rugby(1,:)))/2,(max(rugby(2,:))+min(rugby(2,:)))/2);
        plot2 = fill(rugby(1,:),rugby(2,:),'r');
        axis([-0.5-0.005*d,d*1.005+0.5,-0.5-0.005*d,d*1.005+0.5])
        drawnow
        
        pause(pause_time/10)
    end
else
    no_kick = 1;
end

%%% Response to final ball position %%%
if (min(rugby(2,:)) <= max(goal(2,:)) && (max(rugby(2,:)) >= min(goal(2,:)))) % The ball hits the goal
    pause(0.8)
    if Practise_mode == 1
        set(menu_button,'Visible','off')
        score = 0;
        uicontrol('Style','pushbutton', ...
            'String', 'Goal! Click to retry', ...
            'Position', [225 215 110 24],...
            'Callback', @Retry_Callback);
        uicontrol('Style','pushbutton', ...
            'String', 'Main Menu', ...
            'Position', [240 185 80 24],...
            'Callback', @Menu_Callback);
        quit_button = uicontrol('Style','pushbutton', ...
            'String', 'Quit', ...
            'Position', [240 155 80 24],...
            'Callback', @Quit_Callback);
    else
        score = score + 1;
        if difficulty == 1
            distance_global = distance_global + 5;
        else
            distance_global = distance_global + 10;
        end
        SpheroidKicker2000
    end
elseif  no_kick == 1 || path(k3) < 0 % The ball hits the ground
    pause(0.5)
    set(menu_button,'Visible','off')
    if Practise_mode == 1
        score = 0;
        uicontrol('Style','pushbutton', ...
            'String', 'Miss! Click to retry', ...
            'Position', [225 215 110 24],...
            'Callback', @Retry_Callback);
        uicontrol('Style','pushbutton', ...
            'String', 'Main Menu', ...
            'Position', [240 185 80 24],...
            'Callback', @Menu_Callback);
        quit_button = uicontrol('Style','pushbutton', ...
            'String', 'Quit', ...
            'Position', [240 155 80 24],...
            'Callback', @Quit_Callback);
    else
        set(Score_display,'Visible','off')
        if difficulty == 2 && score > 0
            annotation('textbox',[.385 .35 .2 .3],...
                'String',sprintf('  Score:  %d  (Hard) ',score),...
                'FitBoxToText','on','Color','k','EdgeColor','w');
        else
            annotation('textbox',[.425 .35 .2 .3],...
                'String',sprintf('  Score:  %d    ',score),...
                'FitBoxToText','on','Color','k','EdgeColor','w');
        end
        uicontrol('Style','pushbutton', ...
            'String', 'Main Menu', ...
            'Position', [240 215 80 24],...
            'Callback', @Menu_Callback);
        quit_button = uicontrol('Style','pushbutton', ...
            'String', 'Quit', ...
            'Position', [240 185 80 24],...
            'Callback', @Quit_Callback);
        score = 0;
    end
else % The ball goes past the goal
    set(plot2,'Visible','off')
    pause(0.5)
    set(menu_button,'Visible','off')
    if Practise_mode == 1
        score = 0;
        uicontrol('Style','pushbutton', ...
            'String', 'Miss! Click to retry', ...
            'Position', [225 215 110 24],...
            'Callback', @Retry_Callback);
        uicontrol('Style','pushbutton', ...
            'String', 'Main Menu', ...
            'Position', [240 185 80 24],...
            'Callback', @Menu_Callback);
        quit_button = uicontrol('Style','pushbutton', ...
            'String', 'Quit', ...
            'Position', [240 155 80 24],...
            'Callback', @Quit_Callback);
    else
        set(Score_display,'Visible','off')
        if difficulty == 2 && score > 0
            annotation('textbox',[.385 .35 .2 .3],...
                'String',sprintf('  Score:  %d  (Hard) ',score),...
                'FitBoxToText','on','Color','k','EdgeColor','w');
        else
            annotation('textbox',[.425 .35 .2 .3],...
                'String',sprintf('  Score:  %d    ',score),...
                'FitBoxToText','on','Color','k','EdgeColor','w');
        end
        uicontrol('Style','pushbutton', ...
            'String', 'Main Menu', ...
            'Position', [240 215 80 24],...
            'Callback', @Menu_Callback);
        quit_button = uicontrol('Style','pushbutton', ...
            'String', 'Quit', ...
            'Position', [240 185 80 24],...
            'Callback', @Quit_Callback);
        score = 0;
    end
end

%%% uicontrol click response %%%
    function key_press_Fcn(~,~)
        key_press = 1;
    end

%%% Shape matrix animator manipulations %%%
    function translated_shape = translate(matrix,x_translation,y_translation)
        translated_shape = matrix+[x_translation;y_translation];
    end

    function newshape = rotate(matrix,angle)
        newshape = [cos(angle) -sin(angle); sin(angle) cos(angle)]*matrix;
    end

    function rotated_shape = rotateabout(matrix,rotation,x_center,y_center)
        rotated_shape = translate(matrix,-x_center,-y_center);
        rotated_shape = rotate(rotated_shape,rotation);
        rotated_shape = translate(rotated_shape,x_center,y_center);
    end
end


