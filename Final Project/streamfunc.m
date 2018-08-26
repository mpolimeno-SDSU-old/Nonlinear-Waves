clc 
clear all
dt = .5;
tmax = 70;
t = 0:dt:tmax;
tlength = length(t);
N = 128;
KT = 2*N;

v = 0.0001;

L = 1;
x = linspace(-L,L,KT+1);
x = x(1:KT);
y = x;
[X,Y] = meshgrid(x,y);

Dd1 = 0:N-1;
Dd1(1,1) = 10^(-6);
Dd2 = -N:-1;
Dds = 2*1i*pi/2*[Dd1 Dd2]';

Dy = kron(Dds,ones(KT,1));
Dx = kron(ones(KT,1),Dds);

Dy2 = Dy.^2;
Dx2 = Dx.^2;

w = w0(KT,x,y,X,Y);
wF = fft2(w);

sf = poisson_attempt3(w,KT);
sfF = fft2(sf);

wF = reshape(wF',KT^2,1);

sfF = reshape(sfF',KT^2,1);

streamfunction = zeros(KT,KT,tlength);
streamfunction(:,:,1) = sf;
vorticity = zeros(KT,KT,tlength);
vorticity(:,:,1) = w;

tic
for i = 1:tlength

%Crank-Nicolson
    Dxs = real(ifft2(reshape((Dx.*sfF).',KT,KT).')); 
    Dyw = real(ifft2(reshape((Dy.*wF).',KT,KT).')); 
    Dys = real(ifft2(reshape((Dy.*sfF).',KT,KT).')); 
    Dxw = real(ifft2(reshape((Dx.*wF).',KT,KT).'));
    
    adv = reshape((fft2((Dxs.*Dyw)-(Dys.*Dxw)))',KT^2,1);
    wF1 = (1./((1/dt) - (1/2)*v*(Dx2+Dy2))).*(((1/dt)+(1/2)*v*(Dx2+Dy2)).*wF - adv);
    wF = wF1;
    
    w = real(ifft2(reshape(wF1.',KT,KT).'));
    sf = poisson_attempt3(w,KT);
    streamfunction(:,:,i) = sf;
    vorticity(:,:,i) = w;
    sfF = fft2(sf);
    sfF = reshape(sfF',KT^2,1);
    
end
% v = VideoWriter('badassvortex.avi');
% open(v);

for j = 1:tlength
    pcolor(X,Y,vorticity(:,:,j))
    colorbar;
    shading flat; colormap('jet');
    title(['Advection-Diffusion, t=' num2str(t(j))])
    M(j) = getframe;
%     writeVideo(v,M(j));
end
toc
% close(v);


