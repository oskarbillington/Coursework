% Input: source structure containing attributes .filename, .url, .attributes,
% .formatSpec and .headerlines. Output: table of data from the source. 
% Solved using websave(), textscan() and table(). Downloads from url if 
% filename can't be opened locally.

function file_data = load_data_set(my_source)

% Import data
filename = my_source.filename;
FID = fopen(filename,'r');

if FID == -1 % If the file cannot be opened locally, download it:
    filename = 'tempfile.txt';
    websave(filename, my_source.url);
    FID = fopen(filename,'r');
    sprintf('File was downloaded from the source url because it was not found locally')
end

dataArray = textscan(FID, my_source.formatSpec, 'Headerlines', my_source.headerlines, 'Delimiter', my_source.delimiter, 'MultipleDelimsAsOne', 1, 'WhiteSpace', '', 'ReturnOnError', false); %'TreatAsEmpty', ' ',
fclose(FID);

% Read file into table
file_data = table(dataArray{1:end}, 'VariableNames', my_source.attributes);

% Delete the file and clear temporary variables
if strcmp(filename,'tempfile.txt')
    delete(filename)
end


end