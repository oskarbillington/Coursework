% Input: validation data and weights from regression. Output: mean absolute
% difference between the weighted features and the labels in the test data.


function mae = compute_mean_abs_error(val_D,w)

% Split features from labels
A = val_D(:,1:(end-1)); 
b = val_D(:,end); 

% Approximated labels
approx = A*w;

% Mean square error
mae = mean(abs(b-approx));


end