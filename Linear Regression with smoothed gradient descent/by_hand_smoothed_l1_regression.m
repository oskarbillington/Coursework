% Input: Training data and lambda. Output: Set of weights & t. Fit by
% iterative regression using softmax(?) and gradient descent. Ensures
% numerical stability using the log-exp-trick to shift axis.



function [w,t] = by_hand_smoothed_l1_regression(train_D, lambda)

% INITIAL SETUP:
A = train_D(:,1:(end-1)); % Features
b = train_D(:,end); % Labels

d = size(A,2); % Number of features/dimension of weights vector
n = size(train_D,1); % Number of samples in training set

%syms(sym('w',[1 d]), 'real');
%w = sym('w',[1 d])'; % Create symbolic w-vector to allow easy gradient calculation-implementation

tau = 0.01; % Make as small as possible while keeping numerical stability

% GRADIENT DESCENT:

% Hyperparameters
w = zeros(d,1); % Initialise weights ( (*) zeros ensure numerical stability on first iteration assuming labels don't cause over/underflow on their own) - Actually, they seem to do for tau<0.1... Instead, set g2(1)=labels as A*w(t=1)=0!
step_size = 1e-1; % Adjust manually later
convergence_criterion = 1e-5; % When to stop iterating
max_iterations = 1e3; % Max iterations before aborting and returning error
t = 1; % Count iterations
f_history = zeros(max_iterations,1); % Store function values
function_change = 1; % Just initialise at any value > convergence criterion
f_0 = (lambda/2)*sum(abs(w)) + mean(abs(A*w-b)); % Initial objective function value


disp = figure; % (Just for progress display-plotting)
%axis([-1,100,1,30])
%axis([-1,100,0,1])
hold on
grid on
%set(gca,'YScale','log')
title(sprintf("Lambda = %.1e, Alpha = %.1e",lambda,convergence_criterion))
xlabel('Iterations, t')
ylabel('(Smoothed) Objective Function Value, f_t_a_u(w)')




while (function_change > convergence_criterion) && (t < max_iterations)
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % ADAPTIVE STEP SIZE:
    %if function_change < 100*convergence_criterion % Every time close within convergence:
     %   step_size = step_size*0.90; % Can put this (0.95) outside, call it theta or something
    %end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % SMOOTHED PIECE-WISE LINEAR APPROXIMATION:
    
    % NUMSTAB: Numerical g-values based on previous iteration
    % (non-differentiable, used as constants only)
    g1 = abs(w);
    g2 = abs(b); % @@@@@@@@@@@@ OR abs(A*w + b) ???
    
    % First part of objective function derivative % @@@ can cancel a "tau" in
    % the equation(s) below...
    ddw_l1_norm_w = tau * (exp((-w-g1)/tau) + exp((w-g1)/tau)).^(-1) .* (-exp((-w-g1)/tau)/tau + exp((w-g1)/tau)/tau);
    
    % Second part of objective function derivative
    ddw_mean_absolute_error = mean(tau * (exp((A*w-b-g2)/tau) + exp((b-A*w-g2)/tau)).^(-1) .* (exp((A*w-b-g2)/tau).*A/tau - exp((b-A*w-g2)/tau).*A/tau), 1)';
    
    % Full gradient with respect to each weight
    grad = (lambda/2)*ddw_l1_norm_w + ddw_mean_absolute_error;
    
    %%%%%%%%%%%%%%%%%%%%
    
    % Update the weights
    w = w - grad*step_size;
    
    % Calculate the actual value of the objective function
    if t == 1
        f_prev = f_0;
    else
        f_prev = f_now;
    end
    f_now = (lambda/2)*sum(abs(w)) + mean(abs(A*w-b));
    
    % Calculate the change in objective function-value
    function_change = abs(f_now - f_prev);
    
    % Store the result
    f_history(t) = f_now;
    
    % Plot along the way to see what's up
    plot(t,f_history(t),'b.')
    drawnow
    
    % Next time step
    t = t + 1;
    
end


% Completion message:
if t == max_iterations
    sprintf("Gradient descent iterations aborted due to slow convergence")
    w = 0;
else
    sprintf("Converging objective value: %.4f, lambda: %.1e, steps: %d",f_now,lambda,t)
end

close(disp)

end