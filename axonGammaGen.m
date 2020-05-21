function [x0,y0,axr,L] = axonGammaGen(N,axF,meanD,shp)
% [x0,y0,axr,L] = axonGammaGen(N,axF,meanD,shp) creates a geometry of
% axons with gamma distributed diameters. 
%
% Axons are placed using the description outlined in Hall MG, Alexander
% DC. Convergence and parameter choice for Monte-Carlo simulations of 
% diffusion MRI. IEEE Trans Med Imaging. 2009 Sep;28(9):1354?64. 
% 
% inputs:
%   N = number of axons
%   axF = axon area fraction (range 0-1)
%   meanD = arithmetic mean of axon diameter (um)
%   shp = shape parameter for gamma distribution
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

axD = gamrnd(shp,meanD/shp,[N,1]);
area = sum(pi*(axD/2).^2)/axF;
L = sqrt(area);

axD = sort(axD,'descend');
x0 = [];
y0 = [];
axD0 = [];
rng;
%close all
for n=1:N
    while true
        xi = rand(1)*L;
        yi = rand(1)*L;
        if all(sqrt((xi-x0).^2+(yi-y0).^2) > axD(n)/2+axD0/2)
            break;
        end
    end
    
    % periodic boundary conditions repeat each axon outside the simulation
    % geometry
    x0 = [x0;xi-L;xi-L;xi-L;xi;   xi;xi;   xi+L;xi+L;xi+L];
    y0 = [y0;yi-L;yi;   yi+L;yi-L;yi;yi+L;yi-L;yi;   yi+L];
    axD0 = [axD0;axD(n);axD(n);axD(n);axD(n);axD(n);axD(n);axD(n);axD(n);axD(n)];
end

axr = axD0/2;
