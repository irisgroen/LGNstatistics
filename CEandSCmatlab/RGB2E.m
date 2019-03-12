function [E, El, Ell]= RGB2E(in)
% convert R,G,B to E, El, Ell
R = double(in(:,:,1));
G = double(in(:,:,2));
B = double(in(:,:,3));
% JMG: slightly different values than in PAMI 2001 paper;
%      simply assuming correctly white balanced camera
E   = (0.3*R + 0.58*G + 0.11*B ) / 255.0;
El  = (0.25*R  + 0.25*G - 0.5*B ) / 255.0;
Ell = (0.5*R - 0.5*G) / 255.0;

% % As in original in PAMI 2001 paper;
% E   = (0.06*R + 0.63*G + 0.27*B ) / 255.0;
% El  = (0.3*R  + 0.04*G - 0.35*B ) / 255.0;
% Ell = (0.34*R - 0.6*G+0.17*B) / 255.0;
