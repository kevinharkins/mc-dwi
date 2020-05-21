% this script performs (perhaps incomplete) testing of the mc2sim function

clearvars
close all
clc

now = datetime('now');
filename = sprintf('TEST_mc2sim-%s.md',...
    datetime(now,'Format','yyyy-MM-dd-HH-mm-ss'));
fid = fopen(filename,'w');

if ~fid
    warning('Could not open file for writing. Writing to std out')
    fid = 1;
end

fprintf(fid,'# Testing `mc2sim.m`\n');
fprintf(fid,'\n');
fprintf(fid,'Run at %s\n',now);
fprintf(fid,'\n');

%% jumps from outside an axon near a border


disp('Jump test from outside an axon')
fprintf(fid,'## Jump test from outside an axon\n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 500 000 spins originating from a single point ' ...
    'outside an axon jump in a random direction. About half of the spins ' ...
    'that interact with the outer boundary should reflect back, and not ' ...
    'reflect off the inner boundary.\n']);

N = 500000;
parms.ri = [1];
parms.ro = parms.ri/.65;
parms.xr = [5];
parms.yr = [5];
parms.Lx = 10;
parms.D = [3 3 3]; % [um/ms]
parms.Tn = [100 100];
parms.M0r = 0.5;

dt = .10;
G = zeros(1,1);
A = zeros(size(G));

y0 = parms.yr(1)*ones(N,1);
x0 = parms.xr(1)*ones(N,1)-parms.ro(1)-.3464;

[sig,ph,x,y] = mc2sim(G,A,dt,x0,y0,parms,ones(N,1),zeros(N,1));

xto = parms.xr-parms.ro(1):.0001:parms.xr+parms.ro(1);
xti = parms.xr-parms.ri(1):.0001:parms.xr+parms.ri(1);
figure
plot(x,y,'.')
hold on
plot(x0,y0,'.r');
plot(xto,sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
plot(xto,-sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,-sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
hold off
axis equal

print('-dpng',fullfile('images','jump1.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','jump1.png') ')\n']);
fprintf(fid,'\n');

%% jumps from inside myelin

disp('Jump test from inside myelin')
fprintf(fid,'## Jump test from inside myelin \n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 100 000 spins within a single point inside ' ...
    'myelin jump in a random direction. None of the spins should reflect ' ...
    'back.\n']);

N = 100000;
parms.ri = [1];
parms.ro = parms.ri/.65;
parms.xr = [5];
parms.yr = [5];
parms.Lx = 10;
parms.D = [3 3 3]; % [um/ms]
parms.Tn = [100 100];
parms.M0r = 0.5;

dt = .1;
G = zeros(1,1);
A = zeros(size(G));

y0 = parms.yr(1)*ones(N,1);
x0 = parms.xr(1)*ones(N,1)-parms.ro(1)+.3;

[sig,ph,x,y] = mc2sim(G,A,dt,x0,y0,parms,ones(N,1),zeros(N,1));

xto = parms.xr-parms.ro(1):.0001:parms.xr+parms.ro(1);
xti = parms.xr-parms.ri(1):.0001:parms.xr+parms.ri(1);
figure
plot(x,y,'.')
hold on
plot(x0,y0,'.r');
plot(xto,sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
plot(xto,-sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,-sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
hold off
axis equal

print('-dpng',fullfile('images','jump2.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','jump2.png') ')\n']);
fprintf(fid,'\n');

%% jumps from inside an axon near a border

disp('Jump test from inside an axon')
fprintf(fid,'## Jump test from inside an axon \n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 100 000 spins within a single point inside ' ...
    'an axon jump in a random direction. About half of the spins that ' ...
    'interact with the inner boundary should reflect back, while the ' ...
    'spins that interact with the outer boundary should not reflect.\n']);

N = 100000;
parms.ri = [1];
parms.ro = parms.ri/.65;
parms.xr = [5];
parms.yr = [5];
parms.Lx = 10;
parms.D = [3 3 3]; % [um/ms]
parms.Tn = [100 100];
parms.M0r = .5;

dt = .1;
G = zeros(1,1);
A = zeros(size(G));

y0 = parms.yr(1)*ones(N,1);
x0 = parms.xr(1)*ones(N,1)-parms.ri(1)+.2;

[sig,ph,x,y] = mc2sim(G,A,dt,x0,y0,parms,ones(N,1),zeros(N,1));

xto = parms.xr-parms.ro(1):.0001:parms.xr+parms.ro(1);
xti = parms.xr-parms.ri(1):.0001:parms.xr+parms.ri(1);
figure
plot(x,y,'.')
hold on
plot(x0,y0,'.r');
plot(xto,sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
plot(xto,-sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,-sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
hold off
axis equal

print('-dpng',fullfile('images','jump3.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','jump3.png') ')\n']);
fprintf(fid,'\n');


%% diffusion restricted to intra-axonal space

disp('Intra-axonal restriction')
fprintf(fid,'## Intra-axonal restriction\n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 500 000 spins diffuse inside an axon ' ...
    'for 100 ms. At the end of the simulation, all the spins should ' ...
    'still be inside the axon.\n']);

N = 500000;
parms.ri = [.5];
parms.ro = parms.ri/.65;
parms.xr = [5];
parms.yr = [5];
parms.Lx = 10;
parms.D = [3 3 3]; % [um/ms]
parms.Tn = [100 100];
parms.M0r = 0;

dt = .01;
G = zeros(10000,1);
A = zeros(size(G));

y0 = parms.yr(1)*ones(N,1);
x0 = parms.xr(1)*ones(N,1);

[sig,ph,x,y] = mc2sim(G,A,dt,x0,y0,parms,ones(N,1),zeros(N,1));

xto = parms.xr-parms.ro(1):.0001:parms.xr+parms.ro(1);
xti = parms.xr-parms.ri(1):.0001:parms.xr+parms.ri(1);
figure
plot(x,y,'.')
hold on
plot(x0,y0,'.r');
plot(xto,sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
plot(xto,-sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,-sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
hold off
axis equal

print('-dpng',fullfile('images','intraaxonal.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','intraaxonal.png') ')\n']);
fprintf(fid,'\n');

%% diffusion restricted to intra-myelinic space

disp('Intra-myelin restriction')
fprintf(fid,'## Intra-myelin restriction\n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 100 000 spins diffuse inside myelin ' ...
    'for 100 ms. At the end of the simulation, all the spins should ' ...
    'still be inside myelin.\n']);

N = 100000;
parms.ri = [1];
parms.ro = parms.ri/.65;
parms.xr = [5];
parms.yr = [5];
parms.Lx = 10;
parms.D = [3 3 3]; % [um/ms]
parms.Tn = [100 100];
parms.M0r = inf;

dt = .01;
G = zeros(10000,1);
A = zeros(size(G));

y0 = parms.yr(1)*ones(N,1);
x0 = parms.xr(1)*ones(N,1)-parms.ro(1)+.2;

[sig,ph,x,y] = mc2sim(G,A,dt,x0,y0,parms,ones(N,1),zeros(N,1));

xto = parms.xr-parms.ro(1):.0001:parms.xr+parms.ro(1);
xti = parms.xr-parms.ri(1):.0001:parms.xr+parms.ri(1);
figure
plot(x,y,'.')
hold on
plot(x0,y0,'.r');
plot(xto,sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
plot(xto,-sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,-sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
hold off
axis equal

print('-dpng',fullfile('images','myelin.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','myelin.png') ')\n']);
fprintf(fid,'\n');

%% diffusion restricted to extra-axonal space

disp('Extra-axonal restriction')
fprintf(fid,'## Extra-axonal restriction\n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 500 000 spins diffuse in extra-axonal space ' ...
    'for 100 ms. At the end of the simulation, all the spins should ' ...
    'still be outside of the axons.\n']);

N = 500000;
[xr,yr,ro,Lx] = axonGammaGen(200,0.7,2/0.7,2.331);
parms.ri = 1;
parms.ro = parms.ri/.65;
parms.xr = 2;
parms.yr = 2;
parms.Lx = 4;
parms.D = [3 3 3]; % [um/ms]
parms.Tn = [100 100];
parms.M0r = 0;

dt = .01;
G = zeros(1000,1);
A = zeros(size(G));

y0 = parms.yr(1)*ones(N,1);
x0 = parms.xr(1)*ones(N,1)-parms.ro(1)-.2;
[x0,y0] = axonSeed_ext(parms,N);

[sig,ph,x,y] = mc2sim(G,A,dt,x0,y0,parms,ones(N,1),zeros(N,1));

xto = parms.xr-parms.ro(1):.0001:parms.xr+parms.ro(1);
xti = parms.xr-parms.ri(1):.0001:parms.xr+parms.ri(1);

figure
clf
plot(x,y,'.')
hold on
plot(xto,sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
plot(xto,-sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,-sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
hold off
axis equal

print('-dpng',fullfile('images','extraxonal.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','extraxonal.png') ')\n']);
fprintf(fid,'\n');

%% anisotropic diffusion within myelin

disp('Anisotropic diffusion in myelin')
fprintf(fid,'## Anisotropic diffusion in myelin\n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 100 000 spins originating from the same point ' ...
    'inside myelin anisotropically jump in a random direction. This should ' ...
    'creates a bent ellipse-like shape. \n']);

N = 100000;
parms.ri = [1];
parms.ro = parms.ri/.65;
parms.xr = [5];
parms.yr = [5];
parms.Lx = 10;
parms.D = [.003 3 3]; % [um/ms]
parms.Tn = [100 100];
parms.M0r = .5;

dt = .1;
G = zeros(1,1);
A = zeros(size(G));

y0 = parms.yr(1)*ones(N,1);
x0 = parms.xr(1)*ones(N,1)-parms.ro(1)+.2;

[sig,ph,x,y] = mc2sim(G,A,dt,x0,y0,parms,ones(N,1),zeros(N,1));

xto = parms.xr-parms.ro(1):.0001:parms.xr+parms.ro(1);
xti = parms.xr-parms.ri(1):.0001:parms.xr+parms.ri(1);
figure
plot(x,y,'.')
hold on
plot(x0,y0,'.r');
plot(xto,sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
plot(xto,-sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,-sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
hold off
axis equal

print('-dpng',fullfile('images','anisotropic.png'));

fprintf(fid,'\n');
fprintf(fid,['![](' fullfile('images','anisotropic.png') ')\n']);
fprintf(fid,'\n');

%% equilibrium of spins between different compartments

disp('Equilibrium of spins in each compartment')
fprintf(fid,'## Equilibrium of spins in each compartment\n');
fprintf(fid,'\n');
fprintf(fid,['In this test, 50 000 spins are seeded throughout a geometry ' ...
    'of axons and diffuse for 1 second. At the end of the simulation, ' ...
    'the proportion of spins in each space should be approximately the same ' ...
    'as the start.\n']);

N = 50000;
parms.ri = [2];
parms.ro = [5];%parms.ri/.05;
parms.xr = [10];
parms.yr = [10];
parms.Lx = 20;
parms.D = [.001 3 3]; % [um/ms]
parms.Tn = [100 100];
parms.M0r = 0.5;

dt = .1;
G = zeros(10000,1);
A = zeros(size(G));

[x0,y0] = axonSeed(parms,parms.M0r,N);

[sig,ph,x,y] = mc2sim(G,A,dt,x0,y0,parms,ones(N,1),zeros(N,1));

xt = parms.xr-parms.ro(1):.0001:parms.xr+parms.ro(1);

xto = parms.xr-parms.ro(1):.0001:parms.xr+parms.ro(1);
xti = parms.xr-parms.ri(1):.0001:parms.xr+parms.ri(1);
figure
plot(x,y,'.')
hold on
plot(xto,sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
plot(xto,-sqrt(parms.ro(1).^2-(xto-parms.xr(1)).^2)+parms.yr(1),'r-',...
    xti,-sqrt(parms.ri(1).^2-(xti-parms.xr(1)).^2)+parms.yr(1),'r-')
hold off
axis equal

iarea = pi*parms.ri(1).^2;
marea = pi*parms.ro(1).^2-pi*parms.ri(1).^2;
earea = parms.Lx.^2-pi*parms.ro(1).^2;

% signal fraction of i,m,e space
isf = iarea/(iarea + parms.M0r*marea + earea);
msf = parms.M0r*marea/(iarea + parms.M0r*marea + earea);
esf = earea/(iarea + parms.M0r*marea + earea);

% ideal spin density in i,m,e space
idi = N*isf/iarea;
mdi = N*msf/marea;
edi = N*esf/earea;

% resulting spin density in i,m,e space
x1 = x-parms.xr(1); y1 = y-parms.yr(1);
idr = sum(sqrt(x1.^2+y1.^2)<parms.ri(1))/iarea;
mdr = sum(sqrt(x1.^2+y1.^2)>parms.ri(1) & sqrt(x1.^2+y1.^2)<parms.ro(1))/marea;
edr = sum(sqrt(x1.^2+y1.^2)>parms.ro(1))/earea;

fprintf(fid,'```\n');
fprintf(fid,'intra-axonal: theoretical density: %.2f, resulting density: %.2f\n', idi,idr);
fprintf(fid,'myelin      : theoretical density: %.2f, resulting density: %.2f\n', mdi,mdr);
fprintf(fid,'extra-axonal: theoretical density: %.2f, resulting density: %.2f\n', edi,edr);
fprintf(fid,'```\n\n');

%%

fclose(fid);
