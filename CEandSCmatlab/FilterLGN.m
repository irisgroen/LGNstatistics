function [E, El, Ell] = FilterLGN(in, sigma)

break_off_sigma = 3.;
filtersize = break_off_sigma*sigma;
x=-filtersize:1:filtersize;

Gauss=1/(sqrt(2 * pi) * sigma)* exp((x.^2)/(-2 * sigma * sigma) );
Gx = (x.^2/sigma^4-1/sigma^2).*Gauss;
Gx = Gx-sum(Gx)/size(x,2);
Gx = Gx/sum(0.5*x.*x.*Gx);

Gy = (x.^2/sigma^4-1/sigma^2).*Gauss;
Gy = Gy-sum(Gy)/size(x,2);
Gy = Gy/sum(0.5*x.*x.*Gy);

if size(in,3) ==1
    
    im = double(in)/double(max(in(:)));
    Ex = conv2padded(im, Gx);
    Ey = conv2padded(im, Gy');
    E = sqrt(Ex.^2+Ey.^2);
    El = [];
    Ell = [];
    
else
    
    [E,El,Ell]=RGB2E(in);
   
    im = E;
    Ex = conv2padded(im, Gx);
    Ey = conv2padded(im, Gy');
    E = sqrt(Ex.^2+Ey.^2);
    
    im = El;
    Elx = conv2padded(im, Gx);
    Ely = conv2padded(im, Gy');
    El = sqrt(Elx.^2+Ely.^2);
    
    im = Ell;
    Ellx = conv2padded(im, Gx);
    Elly= conv2padded(im, Gy');
    Ell = sqrt(Ellx.^2+Elly.^2);
end

