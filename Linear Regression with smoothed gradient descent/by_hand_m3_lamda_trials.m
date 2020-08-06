% Gradient descent for smoothed functions
tic
format long


% 1. LOAD DATASET
Source.filename = 'CaliforniaHousingDataSet.txt';
Source.url = 'http://mpawankumar.info/teaching/b1/CaliforniaHousingDataSet.txt'; % Dataset source url
Source.attributes = {'F1','F2','F3','F4','F5','F6','F7','F8','Label'}; % Column descriptions
Source.formatSpec = repmat('%f',[1 9]); % Describe how to read each row (help textscan)
Source.delimiter = '\t'; % Leave as '' if not relevant
Source.headerlines = 0; % Number of lines to skip at the top of the document

Data = load_data_set(Source); % Load the dataset from the source url into a 2D table
Data = table2array(Data); 

%Data = Data(1:100,:); % For dev only: shrinking the dataset


% 4. GRADIENT DESCENT FOR SMOOTHED FUNCTIONS

% Training, validation, and test data sets (1)
frac1 = 0.75;
[trainval_D, test_D] = random_split(Data, frac1);
[test_D, Data_means0, Data_stds0] = normalise_data(test_D, 0);


% Use regression to obtain weights for each lambda trial value
lambda_trials = {0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0};
nt = length(lambda_trials);
n = 5; % Average over 5 iterations per lambda value
history = zeros(n*nt,3);

for x = 1:n
    for i = 1:nt
        
        % Training, validation, and test data sets (2)
        frac2 = 0.8;
        [train_D, val_D] = random_split(trainval_D, frac2);
        [train_D, Data_means1, Data_stds1] = normalise_data(train_D, 0);
        [val_D, Data_means2, Data_stds2] = normalise_data(val_D, 0);

        % Calculate weights
        [weights,iterations] = by_hand_smoothed_l1_regression(train_D, lambda_trials{i});
        
        % Evaluate performance
        mae = compute_mean_abs_error(val_D,weights);
        
        % History, performance logging
        history(i+(x-1)*nt,:) = [lambda_trials{i}, mae, iterations];
        
    end
end


% Use best lambda to get weights over trainval_D and finally test on test_D

% PROCESSING: Plots
%figure
%plot(test_D(:,end),test_D(:,1:(end-1))*weights, 'x', 'Color', 'blue')
%legend("Predictions", "Location", "northwest")
%title('Smoothed Gradient Descent')
%xlabel('Labels')
%hold on
%fplot(@(x) x)
%axis equal
%drawnow





