function [sigN,phN,xN,yN] = mc2sim(G,A,dt,x0,y0,parms,sig,ph)
% [sigN,phN,xN,yN] = mc2sim(G,A,dt,x0,y0,parms,sig,ph)
% 
% Description:
%   Wrapper for CUDA based monte carlo simulation of diffusion-weighted
%   signal
%
% inputs:
%   G - gradient waveform as a 1D vector (mT/m)
%   A - integral of gradient waveform as a 1D vector (mT*ms/m)
%   dt - timestep (ms)
%   x0 - vector of the initial x-component of seed locations (um)
%   y0 - vector of the initial y-component of seed locations (um)
%   parms - structure of tissue parameters
%       parms.ri - inner axon radius (um)
%       parms.ro - outer axon radius (um)
%       parms.xr - center of axon in X (um)
%       parms.yr - center of axon in Y (um)
%       parms.Lx - length of the square geometry (um)
%       parms.D - [Dmr Dmc D0], diffusion characteristics of myelin
%           (radial and circumfrential) and free water (um^2/ms)
%       parms.T1 - [T1m T1ie], T1 of myelin and non-myelin compartments (ms)
%       parms.T2 - [T2m T2ie], T2 of myelin and non-myelin compartments (ms)
%       parms.M0r - ratio of myelin spin density to non-myelin 
%   sig - initial magnitude component of spins (corresponding to x0,y0)
%   ph - initial phase component of spins (corresponding to x0,y0, and sig)
%
% outputs:
%   sigN - vector of the resulting magnitude component of spins
%   phN - vector of the resulting phase component of spins
%   xN - vector of the resulting x-component of seed locations
%   yN - vector of the resulting y-component of seed locations
%
% by Kevin Harkins (kevin.harkins@vanderbilt.edu)


% specified parameters
N = length(x0);
ri = parms.ri; % inner radius
ro = parms.ro; % outer radius
xr = parms.xr;
yr = parms.yr;
Lx = parms.Lx;
D = parms.D; % [um/ms]
R2m = dt/parms.Tn(1);
R2o = dt/parms.Tn(2);
M0r = parms.M0r; % M0 ratio
sd = sqrt(4*dt*D);

% remove axons that are outside of the arena edge + 1 time step
axMask = xr+ro>-max(sd) & xr-ro<Lx+max(sd) & yr+ro>-max(sd) & yr-ro<Lx+max(sd);
ri = ri(axMask);
ro = ro(axMask);
xr = xr(axMask);
yr = yr(axMask);

% transfer data to the gpu
xi = gpuArray(single(x0));
yi = gpuArray(single(y0));
sigi = gpuArray(single(sig));
phi = gpuArray(single(ph));

% setup cuda code, using the default gpu device
BlockMax = 512;
fun = parallel.gpu.CUDAKernel(fullfile('cu','mc2.ptx'), ...
    fullfile('cu','mc2.cu'));
fun.ThreadBlockSize = [min(N,BlockMax) 1 1];
fun.GridSize = [ceil(N/BlockMax) 1];
fun.SharedMemorySize = 0;
g = gpuDevice();

% initalize the random number generator
parallel.gpu.rng('shuffle');

% run simulation
fprintf('%3.0f%% done',0);
x = [xi(:);yi(:)];
sig = [sigi(:);phi(:)];
warning('off','parallel:gpu:kernel:NullPointer');
for p=1:length(G)
    ra = gpuArray.rand(2*N,1,'single');
    wait(g)
    [x, sig] = feval(fun, x, sig, ra, sd(1), sd(2), sd(3), ...
        G(p)*dt*267.5/10000, A(p)*267.7/10000, [ro(:);ri(:)], ...
        [xr(:);yr(:)], length(xr), Lx, M0r, exp(-R2m), exp(-R2o),N);
    wait(g)
    fprintf('\b\b\b\b\b\b\b\b\b%3.0f%% done',p/length(G)*100)
end
fprintf('\n');
warning('on','parallel:gpu:kernel:NullPointer');

% transfer data back to the cpu
xN = gather(x(1:end/2));
yN = gather(x(end/2+1:end));
sigN = gather(sig(1:end/2));
phN = gather(sig(end/2+1:end));

