function [xi,yi] = axonSeed(dist,M0r,N)
% [xi,yi] = axonSeed(dist,M0r,N)
% 
% Description:
%   Wrapper for CUDA code to seed a set of spins onto a distribution of 
%   myelinated axons
%
% inputs:
%   dist - structure of axon distribution
%       dist.ri - inner axon radius
%       dist.ro - outer axon radius
%       dist.xr - center of axon in X
%       dist.yr - center of axon in Y
%       dist.Lx - length of geometry (in X & Y)
%   M0r - ratio of myelin spin density to non-myelin
%   N - number of spins to seed
%
% outputs:
%   xi - Nx1 vector of the x component of seed locations
%   yi - Nx1 vector of the y component of seed locations
%
% by Kevin Harkins (kevin.harkins@vanderbilt.edu)

parallel.gpu.rng('shuffle');
r1 = gpuArray.rand(2*N,1,'single');
r2 = gpuArray.rand(N,1,'single');
x1 = gpuArray.zeros(2*N,1,'single');

BlockMax = 256;
fun = parallel.gpu.CUDAKernel(fullfile('cu','mcSeed.ptx'), ...
    fullfile('cu','mcSeed.cu'));
fun.ThreadBlockSize = [min(N,BlockMax) 1 1];
fun.GridSize = [ceil(N/BlockMax) 1];
fun.SharedMemorySize = 0;

warning('off','parallel:gpu:kernel:NullPointer');

[x] = feval(fun, x1, r1, r2, [dist.ro(:);dist.ri(:)], ...
    [dist.xr(:);dist.yr(:)], length(dist.yr(:)), dist.Lx, M0r, N);

warning('on','parallel:gpu:kernel:NullPointer');

x2 = gather(x);
xi = x2(1:end/2);
yi = x2(end/2+1:end);