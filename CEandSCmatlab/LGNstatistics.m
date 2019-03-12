function [CE, SC, Beta, Gamma] = LGNstatistics(im, viewingdist, dotpitch, fovbeta, fovgamma)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION
% This function calculates the local contrast distribution of an input
% image. Local contrast is computed using multiscale filters that are
% intended to mimic LGN receptive fields. The local contrast magnitude map
% across the image is summarized in a histogram, which is then
% characterized in to ways: 1) using on a Weibull fit (Beta and Gamma
% parameters) and 2) using the mean and coefficient of variation of the
% histogram (contrast energy, CE, and spatial coherence, SC). CE and SC are
% approximations of Beta and Gamma, and are thought to constitute a more
% biologically plausible computation of a local contrast statistic based on
% LGN outputs that the Weibull fit.
%
% INPUT
% im:                     input image (1d (gray) or 3d (color))
% viewingdist:            viewing distance to image (default = 1 meter)
% dotpitch:               monitor dot pitch in meters
% fovBeta:                field of view for estimation of Weibull
%                         beta parameter (default = 1.5 degrees)
% fovGamma:               field of view for estimation of Weibull
%                         gamma parameter (default = 5 degrees)
%
%
% OUTPUT
% CE:                     Contrast energy parameter for 3 color
%                         components (gray , blue-yellow, red-green)
% SC:                     Spatial coherence parameter for 3 color
%                         components (gray , blue-yellow, red-green)
% Beta:                   Weibull beta parameter for 3 color
%                         components (gray , blue-yellow, red-green)
% Gamma:                  Weibull gamma parameter for 3 color
%                         components (gray , blue-yellow, red-green)
%
% ACKNOWLEDGEMENTS
% This code was written by Sennay Ghebreab and modified by H. Steven
% Scholte and Iris Groen. This version of the model was first used in the
% following publication:
%
% Groen IIA, Ghebreab S, Prins H, Lamme VAF, Scholte HS (2013) From image
% statistics to scene gist: Evoked neural activity reveals transition from
% low-level natural image structure to scene category. The Journal of
% Neuroscience 33(48):18813-18824.
%
% Prior publications using a single-scale version of this model are the
% following:
%
% Scholte HS, Ghebreab S, Waldorp L, Smeulders AWM, Lamme VAF (2009)
% Brain responses strongly correlate with Weibull image statistics when
% processing natural images. J Vis 9:1?15. 
%
% Ghebreab S, Smeulders AWM, Scholte HS, Lamme VAF (2009) A biologi-
% cally plausible model for rapid natural image identification. Advances
% in Neural Information Processing Systems 22, pp 629?637.
%
% Please cite these publications if you use this model for your own work.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load ThresholdLGN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default parameters if not specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if  nargin < 5 || isempty(fovgamma)
    fovgamma = 5;
end

if nargin < 4 || isempty(fovbeta)
    fovbeta = 1.5;
end

if nargin < 3 || isempty(dotpitch)
    dotpitch = .35*(10^-3);
end

if nargin < 2 || isempty(viewingdist)
    viewingdist = 1;
end

if nargin < 1
    error('you must specify an input image');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See whether we need to load the image or that we are reading an array %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(im)
    im = imread(im);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set image parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if numel(size(im)) == 2
    IMTYPE = 1; %Gray
elseif numel(size(im)) == 3
    IMTYPE = 2; %Color
end

imsize = [size(im,1) size(im,2)];
minmag1 = double(zeros(imsize));
minmag2 = double(zeros(imsize));
minmag3 = double(zeros(imsize));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters for field of view
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CT0 = 1/75;                         % constant from Geisler&Perry
alpha = (0.106)*1;                  % constant from Geisler&Perry
epsilon2 = 2.3;                     % constant from Geisler&Perry
fovx = round(imsize(2)/2);          % x-pixel loc. of fovea center
fovy = round(imsize(1)/2);          % y-pixel loc. of fovea center
% ex and ey are the x- and y- offsets of each pixel compared to
% the point of focus (fovx,fovy) in pixels.
[ex, ey] = meshgrid(-fovx+1:imsize(2)-fovx,-fovy+1:imsize(1)-fovy);
% eradius is the radial distance between each point and the point
% of gaze.  This is in meters.
eradius = dotpitch .* sqrt(ex.^2+ey.^2);
clear ex ey;
% calculate ec, the eccentricity from the foveal center, for each
% point in the image.  ec is in degrees.
ec = 180*atan(eradius ./ viewingdist)/pi;
% select the pixels that fall within the input visual field of view
imfovbeta = find(ec<fovbeta);
imfovgamma = find(ec<fovgamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Par1 = double(zeros(imsize));
Par2 = double(zeros(imsize));
Par3 = double(zeros(imsize));
Mag1 = double(zeros(imsize));
Mag2 = double(zeros(imsize));
Mag3 = double(zeros(imsize));

offset = 0;
filtnr = 0;
scalenr = 0;
for sigma = [48 24 12 6 3]
    filtnr = filtnr+1;
    
    sigmas = [1 2 4 8 16 32 ];
    v1 = squeeze(ThresholdLGN(:,1));
    t1 = interp1(sigmas,  v1, sigma,'linear');
    v2 = squeeze(ThresholdLGN(:,2));
    t2 = interp1(sigmas, v2, sigma, 'linear');
    v3 = squeeze(ThresholdLGN(:,3));
    t3 = interp1(sigmas, v3, sigma, 'linear');
    
    
    [O1, O2, O3] = FilterLGN(im,sigma);
    [S1] = LocalCOV(O1,sigma);
    E1 = (O1.* max(O1(:)))./(O1+(max(O1(:)).*S1));
    minm1 = E1 - t1;
    index1 = find(minm1 > 0.0000001);
    Par1(index1) = minm1(index1);
    
    
    if IMTYPE == 2
        
        [S2] = LocalCOV(O2,sigma);
        E2 = (O2.* max(O2(:)))./(O2+(max(O2(:)).*S2));
        minm2 = E2 - t2;
        index2 = find(minm2 > 0.0000001);
        Par2(index2) = minm2(index2);
        
        [S3] = LocalCOV(O3,sigma);
        E3 = (O3.* max(O3(:)))./(O3+(max(O3(:)).*S3));
        minm3 = E3 - t3;
        index3 = find(minm3 > 0.0000001);
        Par3(index3) = minm3(index3);        
    end
end

offset = 0;
filtnr = 0;
scalenr = 0;
for sigma = [64 32 16 8 4]
    filtnr = filtnr+1;
   
    
    sigmas = [1 2 4 8 16 32];
    v1 = squeeze(ThresholdLGN(:,1));
    t1 = interp1(sigmas,  v1, sigma,'linear');
    v2 = squeeze(ThresholdLGN(:,2));
    t2 = interp1(sigmas, v2, sigma, 'linear');
    v3 = squeeze(ThresholdLGN(:,3));
    t3 = interp1(sigmas, v3, sigma, 'linear');
    
    
    [O1, O2, O3] = FilterLGN(im,sigma);
    [S1] = LocalCOV(O1,sigma);
    E1 = (O1.* max(O1(:)))./(O1+(max(O1(:)).*S1));
    minm1 = E1 - t1;
    index1 = find(minm1 > 0.0000001);
    Mag1(index1) = minm1(index1);
    
    
    if IMTYPE == 2
        
        [S2] = LocalCOV(O2,sigma);
        E2 = (O2.* max(O2(:)))./(O2+(max(O2(:)).*S2));
        minm1 = E2 - t2;
        index2 = find(minm2 > 0.0000001);
        Mag2(index2) = minm2(index2);
        
        [S3] = LocalCOV(O3,sigma);
        E3 = (O3.* max(O3(:)))./(O3+(max(O3(:)).*S3));
        minm3 = E3 - t3;
        index3 = find(minm3 > 0.0000001);
        Mag3(index3) = minm3(index3);        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Feature Energy and Spatial Coherence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
magnitude = abs(Par1(imfovbeta));
CE(1)= mean(magnitude(:));
if IMTYPE == 2
    magnitude = abs(Par2(imfovbeta));
    CE(2)= mean(magnitude(:));
    magnitude = abs(Par3(imfovbeta));
    CE(3)= mean(magnitude(:));
end

magnitude = abs(Mag1(imfovgamma));
SC(1)= mean(magnitude(:))./std(magnitude(:));
if IMTYPE ==  2
    magnitude = abs(Mag2(imfovgamma));
    SC(2)= mean(magnitude(:))./std(magnitude(:));
    magnitude = abs(Mag3(imfovgamma));
    SC(3)= mean(magnitude(:))./std(magnitude(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Weibull parameters beta and gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nBins = 1000;
magnitude = abs(Par1(imfovbeta));
[ax, h] = createHist(magnitude(:), nBins);
[Beta(1), ~] = weibullMleHist(ax, h);
if IMTYPE == 2
    magnitude = abs(Par2(imfovbeta));
    [ax, h] = createHist(magnitude(:), nBins);
    [Beta(2),~] = weibullMleHist(ax, h);
    magnitude = abs(Par3(imfovbeta));
    [ax, h] = createHist(magnitude(:), nBins);
    [Beta(3),~] = weibullMleHist(ax, h);
end


magnitude = abs(Mag1(imfovgamma));
[ax, h] = createHist(magnitude(:), nBins);
[~, Gamma(1)] = weibullMleHist(ax, h);
if IMTYPE == 2
    magnitude = abs(Mag2(imfovgamma));
    [ax, h] = createHist(magnitude(:), nBins);
    [~, Gamma(2)] = weibullMleHist(ax, h);
    magnitude = abs(Mag3(imfovgamma));
    [ax, h] = createHist(magnitude(:), nBins);
    [~, Gamma(3)] = weibullMleHist(ax, h);
end


