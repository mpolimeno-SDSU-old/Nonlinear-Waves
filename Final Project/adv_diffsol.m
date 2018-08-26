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
x = linspace(-L,L,KT);
x = x(1:KT);
y = x;

% x = x(1:KT); %periodic bc
% y = y(1:KT);
[X,Y] = meshgrid(x,y);



Dd1 = 0:N-1;
Dd2 = -N:-1;
Dds = 2*pi/2*[Dd1 Dd2];
Dds(1) = 10^(-6);

[KX,KY] = meshgrid(Dds,Dds);
K = KX.^2 + KY.^2;
size(K)
K2 = reshape(K,KT^2,1);

% [D,x] = cheb(N); y = x; %fixed bc
% [xx,yy] = meshgrid(x(1:N),y(1:N));
% xx = xx(:); yy = yy(:);

Dy = kron(Dds,ones(KT,1));
Dx = kron(ones(KT,1),Dds);

Dy2 = Dy.^2;
Dx2 = Dx.^2;

w = w0(KT,x,y,X,Y);
wF = fft2(w);

% w = p_fbc(wi,N);

% sf = poisson_attempt3(w,KT);
% sfF = fft2(sf);

wF = reshape(wF,KT^2,1);

[t1, usol] = ode45('avd_diff',t,wF,[],v,K,K2,KT,KX,KY);
v = VideoWriter('vorticity.avi');
open(v);

for i = 1:tlength
    usolr = real(ifft2(reshape(usol(i,:),KT,KT)));
    pcolor(X,Y,usolr)
    colorbar;
    shading flat; colormap('jet');
    title(['Advection-Diffusion, t=' num2str(t(i))])
    M(i) = getframe;
    writeVideo(v,M(i));
end

close(v);

% for j = 1:tlength
%     usoli = real(ifft2(reshape(usol())))
%     pcolor(X,Y,vorticity(:,:,j))
%     colorbar;
%     shading flat; colormap('jet');
%     title(['Advection-Diffusion, t=' num2str(t(j))])
%     M(j) = getframe;
% %     writeVideo(v,M(j));
% end


% sfF = reshape(sfF',KT^2,1);
% 
% streamfunction = zeros(KT,KT,tlength);
% streamfunction(:,:,1) = sf;
% vorticity = zeros(KT,KT,tlength);
% vorticity(:,:,1) = w;
% 
% lop = (ones(KT^2,1)./(ones(KT^2,1)-dt*v*(Dx2+Dy2)));
% 
% dxsf = real(ifft2(reshape((Dx.*sfF).',KT,KT).'));
% dyw = real(ifft2(reshape((Dy.*wF).',KT,KT).'));
% dysf = real(ifft2(reshape((Dy.*sfF).',KT,KT).'));
% dxw = real(ifft2(reshape((Dx.*wF).',KT,KT).'));
% adv = dxsf.*dyw-dysf.*dxw;
% adv = fft2(reshape(adv',KT^2,1));
% wFm1 = wF;
% wFm2 = wFm1;
% wFm3 = wFm2;
% 
% % advm1 = adv;
% % advm2 = advm1;
% % advm3 = advm2;
% % tstep = zeros(length(Dx),1);
% for i = 1:tlength
%     %     adam-bashford 4 
% %     dxsf = real(ifft2(reshape((Dx.*sfF).',KT,KT).'));
% %     dyw = real(ifft2(reshape((Dy.*wF).',KT,KT).'));
% %     dysf = real(ifft2(reshape((Dy.*sfF).',KT,KT).'));
% %     dxw = real(ifft2(reshape((Dx.*wF).',KT,KT).'));
% %     adv = dxsf.*dyw-dysf.*dxw;
% %     adv = fft2(reshape(adv',KT^2,1));
% %     wF1 = lop.*(wF + wFm1/3 + dt*((55/24)*adv-(59/24)*(advm1)...
% %     +(37/24)*(advm2)-(3/8)*advm3));
% %     wFm1 = wF;
% %     wF = wF1;
% %     
% %     advm3 = advm2;
% %     advm2 = advm1;
% %     advm1 = adv;
% %     adam-moulton 3
% %     lop1 = ifft2(reshape((Dy.*sfF.*Dx).',KT,KT).');
% %     lop2 = ifft2(reshape((Dx.*sfF.*Dy).',KT,KT).');
% %     lop3 = reshape(fft2(lop1-lop2)',KT^2,1);
% %     
% %     lop = v*(Dx2+Dy2)-(lop3);
% %     lop = v*(Dx2+Dy2)-(-(Dy.*sfF.*Dx)+(Dx.*sfF.*Dy));
% %     invlop = ones(KT^2,1)./(ones(KT^2,1)-dt*(3/8)*lop);
% %     wF1 = invlop.*(wF + dt*lop.*((19/24)*wFm1-(5/24)*wFm2+(1/24)*wFm3));
% %     wFm1 = wF;
% %     wF = wF1;
% %     
% %     wFm3 = wFm2;
% %     wFm2 = wFm1;
% %     wFm1 = wF;
% 
%     wF1 = avd_diff(0,wF,KT,sfF,v,Dx,Dy,Dx2,Dy2);
%     
%     
% 
%     
% %     wF = wF1;
% 
% %     [t1,wF1] = ode45(@avd_diff, [0 dt], wF,[],KT,sfF,v,Dx,Dy,Dx2,Dy2); %KT-4 
% 
% %Crank-Nicolson
% %     Dxs = real(ifft2(reshape((Dx.*sfF).',KT,KT).')); 
% %     Dyw = real(ifft2(reshape((Dy.*wF).',KT,KT).')); 
% %     Dys = real(ifft2(reshape((Dy.*sfF).',KT,KT).')); 
% %     Dxw = real(ifft2(reshape((Dx.*wF).',KT,KT).'));
% %     
% %     adv = reshape((fft2((Dys.*Dxw)-(Dxs.*Dyw)))',KT^2,1);
% %     
% %     wF1 = (1./((1/dt) - (1/2)*v*(Dx2+Dy2))).*(((1/dt)+(1/2)*v*(Dx2+Dy2)).*wF - adv);
% 
% %     wF1 = wF + dt*(avd_diff(0,wF,KT,sfF,v,Dx,Dy,Dx2,Dy2));
% %     
% %     w = real(ifft2(reshape(wF1.',KT,KT).'));
% % %     
% %     wF = wF1;
% %     wF = fft2(w);
% %     wF = wF(:);
%     
%     w = real(ifft2(reshape(wF1.',KT,KT).'));
% %     wF = wF1(end,:);
%     sf = poisson_attempt3(w,KT);
%     streamfunction(:,:,i) = sf;
%     vorticity(:,:,i) = w;
%     sfF = fft2(sf);
%     sfF = reshape(sfF',KT^2,1);
%     
% end
% 
% % v = VideoWriter('badassvortex.avi');
% % open(v);
% 
% for j = 1:tlength
%     pcolor(X,Y,vorticity(:,:,j))
%     colorbar;
%     shading flat; colormap('jet');
%     title(['Advection-Diffusion, t=' num2str(t(j))])
%     M(j) = getframe;
% %     writeVideo(v,M(j));
% end
% 
% % close(v);


