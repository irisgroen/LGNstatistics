function [CE, SC, Beta, Gamma, filenames] = run_LGNstatistics()

% wrapper function that executes LGNstatistics on each image in the
% directory it is run from

names = dir('*.jpg');
names2 = dir('*.jpeg');
names3 = dir('*.bmp');
names4 = dir('*.png');
names5 = dir('*.tif');
    
names = [names; names2; names3; names4; names5];

CE = nan(length(names),3); 
SC = nan(length(names),3);
Beta = nan(length(names),3);
Gamma = nan(length(names),3);
filenames = cell(size(names)); 

for cNames = 1:length(names)
    if names(cNames).bytes > 0
        im = imread(names(cNames).name);
        [CE(cNames,:),SC(cNames,:), Beta(cNames,:), Gamma(cNames,:)] = LGNstatistics(im);
    end
    filenames{cNames} = names(cNames).name;
    disp(cNames);
end

matname = 'LGNstatistics';

try 
	save(matname,'CE','SC', 'Beta', 'Gamma', 'filenames');
catch
    delete(matname);
	save(matname,'CE','SC', 'Beta', 'Gamma', 'filenames');
end