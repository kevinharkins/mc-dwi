function b = GxtoB(Gx,dt)
% b = GxtoB(Gx,dt) calculates the b-value of a diffusion waveform
%      
%   inputs:
%      Gx = gradient waveform (mT/m)
%      dt = time step (ms)
%      
%   output:
%      b = b-value for diffusion gradient waveform Gx 
%          (assuming Gx integrates to 0) (ms/um^2)
%
% by Kevin Harkins (kevin.harkins@vanderbilt.edu)

Gx_um = Gx/10e3; % mT/um
gamma = 267.513; % rad/ms/mT
b = sum(cumsum(gamma*Gx_um).^2)*dt.^3; % ms/um^2