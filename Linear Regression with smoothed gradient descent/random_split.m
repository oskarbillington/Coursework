% Input: dataset in matrix/table form, assuming each row is one sample where the
% last element is its output/label; and a scalar value representing the
% fraction of the samples returned as training data. Output: two matrices
% containing training data and test data respectively, randomly split and ordered.
% Solved using randperm(). 


function [train_D,test_D] = random_split(D, frac)

% Convert table to array
if isa(D, 'table')
    D = table2array(D);
end

n_samples = size(D,1); % Total number of samples to be distributed
n_train = round(n_samples*frac); % Number of training samples
train_D = zeros(n_train, size(D,2)); % Initialise training data matrix

selection = randperm(n_samples, n_train); % Selecting which samples to train on AND in which order

for k = 1:n_train % Create subsequent training samples:
    train_D(k,:) = D(selection(k),:);
    D(selection(k),1) = " "; % Mark the sample by making first attribute "NaN" but don't remove it yet as it would mess up the rest of this loop
end

test_D = D(~isnan(D(:,1)),:); % Create test data by simply removing all marked rows from D

end

% DISCUSSION:
% Even though selecting out test samples rather than training samples would 
% most likely result in faster computation (because, most likely, frac > 0.5), 
% the following method randomises the order of the training samples, which 
% can impact the optimisation algorithm, whereas leaving the test samples 
% in the original order makes no difference. Making a second loop to
% randomise the samples would be even less efficient.