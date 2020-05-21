function [x0,y0,axr,L] = axonLogNormGen(N,axF,meanD,sigma)
% [x0,y0,axr,L] = axonLogNormGen(N,axF,meanD,sigma) creates a geometry of
% axons with lognormal distributed diameters. 
%
% Axons are placed using the description outlined in Hall MG, Alexander
% DC. Convergence and parameter choice for Monte-Carlo simulations of 
% diffusion MRI. IEEE Trans Med Imaging. 2009 Sep;28(9):1354?64. 
% 
% inputs:
%   N = number of axons
%   axF = axon area fraction (range 0-1)
%   meanD = arithmetic mean of axon diameter (um)
%   sigma = width of the log spaced diameters
%
% outputs: 
%   x0 = vector containing center position of axon in x (um)
%   y0 = vector containing center position of axon in y (um)
%   axr = vector containting outer axon radius (um)
%   L = boundry of simulation (um)
%
% Note: If axon packing is too dense (i.e. axF is too high), this procedure
% will loop inifinitely. 
%
% by Kevin Harkins (kevin.harkins@vanderbilt.edu)

mu = log(meanD)-sigma.^2/2;
axD = lognrnd(mu,sigma,N,1);

area = sum(pi*(axD/2).^2)/axF;
L = sqrt(area);

axD = sort(axD,'descend');
rng;
niter = 0;
maxiter = 100000;
while true
    x0 = [];
    y0 = [];
    axD0 = [];
    for n=1:N
        niter = 0;
        while niter<maxiter
            xi = rand(1)*L;
            yi = rand(1)*L;
            if all(sqrt((xi-x0).^2+(yi-y0).^2) > axD(n)/2+axD0/2)
                break;
            end
            niter = niter+1;
        end
        if niter==maxiter
            break;
        end
        x0 = [x0;xi-L;xi-L;xi-L;xi;   xi;xi;   xi+L;xi+L;xi+L];
        y0 = [y0;yi-L;yi;   yi+L;yi-L;yi;yi+L;yi-L;yi;   yi+L];
        axD0 = [axD0;axD(n);axD(n);axD(n);axD(n);axD(n);axD(n);axD(n);axD(n);axD(n)];
        n;
    end
    if niter==maxiter
        % reshuffle, and start over
        continue;
    end
    break;
end

axr = axD0/2;
