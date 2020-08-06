% Input: data to normalise along columns, & bias (boolean). Output:
% z-normalised data (features and labels), mean and std to enable 
% reversing. Solved using MATLAB's mean() and std().

function [ND,means,stds] = normalise_data(D, bias)

dim = size(D,2);
ND = ones(size(D));
means = zeros(1,dim-bias);
stds = means;

for i = (1+bias):dim % Skip first column if bias == 1
    column = D(:,i); % Grab the column
    s = std(column);
    m = mean(column);
    if s ~= 0 % If all column elements are identical, set them to zero (i.e. subtract mean but don't divide by s)
        z = (column - m)./s; % Normalise
        ND(:,i) = z;
    else
        ND(:,i) = 0;
    end
    means(i-bias) = m;
    stds(i-bias) = s;
end

end