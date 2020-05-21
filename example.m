clearvars
close all
clc

%% define experimental parameters

Nspin = 100000;

te = 100; % experimental echo time, ms
dt = 0.01; % time step, ms

bval = linspace(0,3,13); % ms/um^2, diffusion b-value 
bigDelta = 50; % ms
litDelta = 10; % ms

[Gdiff,t] = makePGwave(dt,te,bigDelta,litDelta);

%% define the simulation arena of axons

Nax = 200; % number of axons to use in the arena
axf = 0.6; % axon volume fraction
axd = 1.0; % um, arithmetic mean axon diameter of the distribution
sigma = 0.5; % standard deviation of the log-spaced diameters
gratio = 0.7; % axonal g-ratio

% create the arena structure
[parms.xr,parms.yr,parms.ro,parms.Lx] = axonLogNormGen(Nax,axf,axd,sigma);
parms.ri = gratio*parms.ro;

%% show the gradient waveform and axon geometry

figure(1)
clf
set(gcf,'PaperPosition',[0 0 6 3])
subplot(121)
plot(t,Gdiff,'linewidth',2)
xlabel('time/ms')
ylabel('normalized gradient strength')
grid on
title('gradient waveform')

subplot(122)
plotDist(parms);
title('simulation arena')

print -dpng images/example_sim.png

%% define tissue parameters

parms.T2 = [15 80]; % T2 in myelin and non-myelin compartments, ms

% [Dmr Dmc Die], the "free" diffusion coefficients in myelin (in the radial 
% and circumferential direction) and non-myelin compartments
parms.D = [0.001 3.0 3.0]; % um^2/ms 

% ratio of the spin density in myelin to intra/extra-axonal space
parms.M0r = 0.5; 

%% run the simulation

parms.dt = dt;
parms.N = Nspin;
[sig,isig,msig,esig] = dwmriWAVE(parms,Gdiff,bval);

%% plot some results

figure
plot(bval,sig,'o-','linewidth',2)
hold on
plot(bval,isig,'o-','linewidth',2)
plot(bval,msig,'o-','linewidth',2)
plot(bval,esig,'o-','linewidth',2)
hold off
grid on

xlabel('b-value/(ms/um^2))')
ylabel('signal')
legend('total','intra-axonal','intra-myelinic','extra-axonal',...
    'location','northeast');

print -dpng images/example_sig.png

