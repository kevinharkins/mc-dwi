% this script performs (perhaps incomplete) testing of the dwmriWAVE function

clearvars
clc
close all

now = datetime('now');
filename = sprintf('TEST_dwmriWAVE-%s.md',...
    datetime(now,'Format','yyyy-MM-dd-HH-mm-ss'));
fid = fopen(filename,'w');

if ~fid
    warning('Could not open file for writing. Writing to std out')
    fid = 1;
end

fprintf(fid,'# Testing `dwmriWAVE.m`\n');
fprintf(fid,'\n');
fprintf(fid,'Run at %s\n',now);
fprintf(fid,'\n');


%% Test T2 relaxation test

fprintf(fid,'## T2 Relaxation\n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 100 000 spins diffuse and relax for up to 50 ms ' ...
    'with T2 = 20 ms. The analytic solution is an exponential decay.\n']);

te = 5:5:40;
parms.T2 = [20 20]; % t2 in myelin and non-myelin compartments
parms.T1 = [1e10 1e10]; % t1 in myelin and non-myelin compartments
parms.D = [3 3 3]; % "free" diffusion coefficients in myelin (radial, circumferential) and non-myelin compartments
parms.M0r = 1.0; % ratio of spin density

parms.N = 100000; % # of spins
parms.dt = 0.02; % timestep

parms.xr = [];
parms.yr = [];
parms.ro = [];
parms.ri = [];
parms.Lx = 100; % um
bvalN = 0;

for n=1:length(te)
    Gdiff = zeros((te(n)/parms.dt),1);
    sig(n) = dwmriWAVE(parms,Gdiff,bvalN);
end

ana = exp(-te./parms.T2(1));

figure
clf
plot(te,sig,'o',te,ana,'-',te,sig-ana,'^')
grid on
xlabel('echo time (ms)')
ylabel('signal')
title('T2 relaxation')
legend('simulated signal','analytic T2 decay','error', ...
    'location','northeast');

print('-dpng',fullfile('images','relax.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','relax.png') ')\n']);
fprintf(fid,'\n');

fprintf(fid,'RMSE: %e\n\n',norm(sig-ana));

%% DWI of free diffusion

fprintf(fid,'## Free diffusion in an empty geometry\n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 100 000 spins diffuse in an empty arena ' ...
    'with bigDelta = 20 ms, littleDelta = 5 ms. Signal vs b-value ' ...
    'should match an exponential decay.\n']);

parms.T2 = [1e10 1e10]; 
parms.T1 = [1e10 1e10]; 
parms.D = [3 3 3]; 
parms.M0r = 1.0; 

parms.N = 100000; 
parms.dt = 0.02; 

parms.xr = [];
parms.yr = [];
parms.ro = [];
parms.ri = [];
parms.Lx = 100; % um

bigD = 20; % ms
litD = 5; % ms
te = 30; % ms
Gdiff = makePGwave(parms.dt,te,bigD,litD);
bvalN = linspace(0,2,11);

sig = dwmriWAVE(parms,Gdiff,bvalN);

ana = exp(-bvalN.*parms.D(3));

figure
plot(bvalN,sig,'o',bvalN,ana,'-',bvalN,sig-ana,'o')
grid on
xlabel('b-value (ms/\mu m^2)')
ylabel('signal')
title('Free diffusion')
legend('simulated signal','analytic diffusion decay','error', ...
    'location','northeast');

print('-dpng',fullfile('images','freediff.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','freediff.png') ')\n']);
fprintf(fid,'\n');

fprintf(fid,'RMSE: %e\n\n',norm(sig-ana));


%% DWI of free diffusion in an unrestricted geometry

fprintf(fid,'## Free diffusion in an unrestricted geometry\n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 100 000 spins diffuse in an arena of axons ' ...
    'that impose no restriction. Signal vs b-value ' ...
    'should still match an exponential decay.\n']);

parms.T2 = [1e10 1e10]; % t2 in myelin and non-myelin compartments
parms.T1 = [1e10 1e10]; % t1 in myelin and non-myelin compartments
parms.D = [3 3 3]; % "free" diffusion coefficients in myelin (radial, circumferential) and non-myelin compartments
parms.M0r = 1.0; % ratio of spin density
parms.N = 100000; % # of spins
parms.dt = 0.02; % timestep

% create a white matter geometry, but it is totally unrestrictive 
Nax = 200;
axf = 0.7;
axd = 2.0;
g = 0.7;
[parms.xr,parms.yr,parms.ro,parms.Lx] = axonLogNormGen(Nax,axf,axd,0.5);
parms.ri = parms.ro*g;

bigD = 20; % ms
litD = 5; % ms
te = 30; % ms
Gdiff = makePGwave(parms.dt,te,bigD,litD);
bvalN = linspace(0,2,11);

sig = dwmriWAVE(parms,Gdiff,bvalN);

ana = exp(-bvalN.*parms.D(3));

figure
plot(bvalN,sig,'o',bvalN,ana,'-',bvalN,sig-ana,'o')
grid on
xlabel('b-value (ms/\mu m^2)')
ylabel('signal')
title('unrestricted diffusion')
legend('simulated signal','analytic diffusion decay','error', ...
    'location','northeast');

print('-dpng',fullfile('images','unresdiff.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','unresdiff.png') ')\n']);
fprintf(fid,'\n');

fprintf(fid,'RMSE: %e\n\n',norm(sig-ana));

figure
plotDist(parms)
title('unrestricted diffusion geometry')

print('-dpng',fullfile('images','unresdiff_arena.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','unresdiff_arena.png') ')\n']);
fprintf(fid,'\n');

%% look at intra-axonal signal 

fprintf(fid,'## Intra-axonal diffusion\n');
fprintf(fid,'\n');
fprintf(fid,['In this test, intra-axonal signal is compared to an ' ...
    'analytic expression given by van Gelderen & DesPres J Magn 1994.\n']);

parms.T2 = [1e10 1e10]; % t2 in myelin and non-myelin compartments
parms.T1 = [1e10 1e10]; % t1 in myelin and non-myelin compartments
parms.D = [1e-10 3 3]; % "free" diffusion coefficients in myelin (radial, circumferential) and non-myelin compartments
parms.M0r = 0.0; % ratio of spin density
parms.N = 100000; % # of spins
parms.dt = 0.01; % timestep

% create a white matter geometry, but it is totally unrestrictive 
Nax = 200;
axf = 0.5;
axd = 12.0;
g = 0.98;
[parms.xr,parms.yr,parms.ro,parms.Lx] = axonDeltaGen(Nax,axf,axd);
parms.ri = parms.ro*g;

bigD = 20; % ms
litD = 5; % ms
te = 30; % ms
Gdiff = makePGwave(parms.dt,te,bigD,litD);
bvalN = linspace(0,2,11);

[sig,intsig] = dwmriWAVE(parms,Gdiff,bvalN);

b = GxtoB(Gdiff,parms.dt);
gmax = sqrt(bvalN/b);
gmax(b==0) = 0;
ana = intsig(1)*vanGelderen(bigD,litD,gmax,axd*g/2,parms.D(3));

figure
plot(bvalN,intsig,'o',bvalN,ana,'-',bvalN,intsig-ana,'^')
grid on
xlabel('b-value (ms/\mu m^2)')
ylabel('intra-axonal signal')

print('-dpng',fullfile('images','intrasig.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','intrasig.png') ')\n']);
fprintf(fid,'\n');
fprintf(fid,'RMSE: %e\n\n',norm(intsig-ana));
fprintf(fid,'Note: the van Gelderen approximation is only valid for low b*Dapp\n');
