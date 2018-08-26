clc 
clear all
dt = .5;
tmax = 50;
t = 0:dt:tmax;
tlength = length(t);
N = 128;
KT = 2*N;

v = 0.0001;

L = 1;
h = L/N;
x = linspace(-L,L,KT+1);
x = x(1:KT);
y = x;

[X,Y] = meshgrid(x,y);

Dd1 = 0:N-1;
Dd1(1,1) = 10^(-6);
Dd2 = -N:-1;
Dds = 2*pi/2*[Dd1 Dd2]';


Dx = kron(Dds,ones(KT,1));
Dy = kron(ones(KT,1),Dds);


Dy2 = Dy.^2;
Dx2 = Dx.^2;

w = w0(KT,x,y,X,Y);
wF = fft2(w);

wF = reshape(wF',KT^2,1);

[t1,wF1] = ode45('adv_diff', t, wF,[],KT,v,Dx,Dy,Dx2,Dy2);

% v = VideoWriter('doublevortex.avi');
% open(v);

for i = 1:tlength
    usolr = real(ifft2(reshape(wF1(i,:).',KT,KT).'));
    pcolor(X,Y,usolr)
    colorbar;
    shading flat; colormap('jet');
    title(['Advection-Diffusion, t=' num2str(t(i))])
    M(i) = getframe;
    writeVideo(v,M(i));
end

% close(v);


