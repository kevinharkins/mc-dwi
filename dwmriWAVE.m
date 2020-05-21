function [dwsig,intsig,myesig,extsig] = dwmriWAVE(parms,Gdiff,bval) 
% [dwsig,intsig,myesig,extsig] = dwmriWAVE(parms,Gdiff,bval)
% 
% Description:
%   Performs a monte carlo simulation of diffusion-weighted signal for a
%   combination of white matter geometry and diffusion-weighted gradients.
%
% inputs:
%   parms - structure of tissue parameters
%       parms.ri - vector of inner axon radius (um)
%       parms.ro - vector of outer axon radius (um)
%       parms.xr - vector of center of axons in X (um)
%       parms.yr - vector of center of axons in Y (um)
%       parms.Lx - scalar length of the square geometry (um) 
%       parms.D - [Dmr Dmc D0], diffusion characteristics of myelin ...
%           (radial and circumfrential) and free water (um^2/ms)
%       parms.T1 - [T1m T1ie], T1 of myelin and non-myelin compartments (ms)
%       parms.T2 - [T2m T2ie], T2 of myelin and non-myelin compartments (ms)
%       parms.M0r - ratio of myelin spin density to non-myelin
%       parms.dt - time step (ms)
%   Gdiff - 2D matrix of diffusion-weighted gradient waveforms, scaled 
%       between [0,1]. The waveforms are automatically scaled to all   
%       combinations of the b-values provided. 
%       size(Gdiff) = [#time steps, # experiments]
%   bval - vector of b-values (ms/um^2)
%
% outputs: (compartments are assigned at the end of the simulation)
%   dwsig - diffusion-weighted signal
%   intsig - intra-axonal signal 
%   myesig - myelin signal
%   extsig - extra-axonal signal
%
% by Kevin Harkins (kevin.harkins@vanderbilt.edu)

dt = parms.dt;

% seed spins
disp('seeding....')
[x0,y0] = axonSeed(parms,parms.M0r,parms.N);
disp('seeding done!')

sig = ones(parms.N,1,'single');
ph = zeros(parms.N,1,'single');

% loop over the number of dw gradient waveforms
fprintf('\n')
for n=1:size(Gdiff,2)
    
    disp(['experiment #' num2str(n)]);
    A1 = cumtrapz(Gdiff(:,n))*dt;
    parms.Tn = parms.T2;
    
    % 
    disp('simulating...')
    [sigN,phN,xN,yN] = mc2sim(Gdiff(:,n),A1,dt,x0,y0,parms,sig,ph);
    
    % mask spins into compartments
    intra = false(size(sigN));
    myelin = false(size(sigN));
    for m=1:length(parms.xr)
        intra = (intra | ...
            (sqrt((xN-parms.xr(m)).^2+(yN-parms.yr(m)).^2) < parms.ri(m)));
        myelin = (myelin | ...
            (sqrt((xN-parms.xr(m)).^2+(yN-parms.yr(m)).^2) < parms.ro(m))) ...
            & ~intra;
    end
    extra = ~(myelin | intra);
    
    % occasionally the GeForce GPUs hiccup, and return invalid values...
    % create a mask for those
    isvalid = ~isnan(phN) & ~isinf(phN) & ~isnan(sigN) & ~isinf(sigN);
    nvalid = sum(isvalid);
    
%     if nvalid<size(sigN,1)
%         fprintf('dropping %i bad values\n', size(sigN,1)-nvalid)
%     end
    
    % calculate gradient amplitudes
    b = GxtoB(Gdiff(:,n),dt);
    if b == 0
        Gmax = 0;
    else
        Gmax = sqrt(bval/b);
    end
    
    % calculate diffusion signals
    dwsig(n,:) = real(sum(bsxfun(@times,gather(sigN(isvalid)),...
        exp(1i*bsxfun(@times,phN(isvalid),Gmax/1))),1)/nvalid);
    intsig(n,:) = real(sum(bsxfun(@times,gather(sigN(intra&isvalid)),...
        exp(1i*bsxfun(@times,phN(intra&isvalid),Gmax/1))),1)/nvalid);
    myesig(n,:) = real(sum(bsxfun(@times,gather(sigN(myelin&isvalid)),...
        exp(1i*bsxfun(@times,phN(myelin&isvalid),Gmax/1))),1)/nvalid);
    extsig(n,:) = real(sum(bsxfun(@times,gather(sigN(extra&isvalid)),...
        exp(1i*bsxfun(@times,phN(extra&isvalid),Gmax/1))),1)/nvalid);
    
    disp(' ');
end
