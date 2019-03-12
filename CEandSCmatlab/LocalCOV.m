function out = LocalCOV(E, sigma)

break_off_sigma = 3.;
filtersize = round(break_off_sigma*sigma);
H = fspecial('gaussian', filtersize, sigma);

term1 = imfilter(E.^2, H, 'replicate');
term2 = imfilter(E, H, 'replicate').^2;
local_std = sqrt(max(term1 - term2, 0));

local_mean = imfilter(E, H, 'replicate')+realmin;

out= local_std./local_mean;

end