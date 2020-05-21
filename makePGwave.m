function [Gdiff,t] = makePGwave(dt,te,bigDelta,litDelta)
% Gdiff = makePGwave(dt,te,bigDelta,litDelta) creates a pulsed gradient
% diffusion waveform with given experimental parameters
% 
% inputs:
%   dt = time step (ms)
%   te = experimental echo time (ms)
%   bigDelta = spacing between diffusion encoding pulses (ms)
%   litDelta = duration of diffusion encoding pulses (ms)
%
% outputs: 
%   Gdiff = envelope of a pulsed gradient diffusion waveform
%
% by Kevin Harkins (kevin.harkins@vanderbilt.edu)

t = (dt:dt:te)';
Gdiff = zeros(size(t));
Gdiff(t>(te/2-bigDelta/2-litDelta/2) & t<(te/2-bigDelta/2+litDelta/2)) = 1;
Gdiff(t>(te/2+bigDelta/2-litDelta/2) & t<=(te/2+bigDelta/2+litDelta/2)) = -1;
