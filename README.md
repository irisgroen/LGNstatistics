LGNstatistics

========= DESCRIPTION ========= 

This is a set of Matlab functions that can be used to compute local contrast statistics for an arbitrary input image. 

Local contrast is computed using multiscale filters that are intended to mimic LGN receptive fields. The local contrast magnitude map across the image is summarized in a histogram, which is then characterized in to ways: 1) using on a Weibull fit (Beta and Gamma parameters) and 2) using the mean and coefficient of variation of the histogram (contrast energy, CE, and spatial coherence, SC). CE and SC are approximations of Beta and Gamma, and are thought to constitute a more biologically plausible computation of a local contrast statistic based on LGN outputs that the Weibull fit.

====== INSTRUCTIONS =======

The main function (LGNstatistics.m) can be called using run_LGNstatistics.m, which is a wrapper function that executes the code for all images within the directory the function is called from. It will take into account image files with the following extensions:

*.jpg
*.jpeg
*.png
*.bmp
*.tif

LGNstatistics.m has default settings describing the experimental setup, but the following settings can be customized:

viewingdist:            viewing distance to image (default = 1 meter)
dotpitch:               monitor dot pitch in meters
fovBeta:                field of view for estimation of CE and Beta (default = 1.5 degrees)
fovGamma:               field of view for estimation of SC and Gamma (default = 5 degrees)

The script generates the following output per image
CE:                     Contrast energy parameter for 3 color
                         components (gray , blue-yellow, red-green)
SC:                     Spatial coherence parameter for 3 color
                         components (gray , blue-yellow, red-green)
Beta:                   Weibull beta parameter for 3 color
                         components (gray , blue-yellow, red-green)
Gamma:                  Weibull gamma parameter for 3 color
                         components (gray , blue-yellow, red-green)

====== ACKNOWLEDGEMENTS =======

This code was written by Sennay Ghebreab and modified by H. Steven Scholte and Iris Groen, at the University of Amsterdam. This version of the model was first used in the following publication:

Groen IIA, Ghebreab S, Prins H, Lamme VAF, Scholte HS (2013) From image statistics to scene gist: Evoked neural activity reveals transition from low-level natural image structure to scene category. The Journal of Neuroscience 33(48):18813-18824.

Prior publications using a single-scale version of this model are the following:

Scholte HS, Ghebreab S, Waldorp L, Smeulders AWM, Lamme VAF (2009) Brain responses strongly correlate with Weibull image statistics when processing natural images. J Vis 9:1-15. 

Ghebreab S, Smeulders AWM, Scholte HS, Lamme VAF (2009) A biologically plausible model for rapid natural image identification. Advances in Neural Information Processing Systems 22, 629-637.

Please cite these publications if you use this model for your own work.

Please contact iris.groen[at]nyu.edu for questions.

