%solve the Poisson Equation
clear all; help fftpoi 

%initialize parameters

eps0 = 8.8542*10^(-12);
N=50;
L=1;
h=L/N;
x=((1:N)-1/2)*h;
y=x;
fprintf('System is a square of length %g \n',L);

%set up charge density

rho= zeros(N,N);
M = input('Enter number of line: '); %this gives you the number of charges that you can input (i.e vorticies)

for i=1:M
    fprintf('\n For charge #%g \n', i);
    r = input('Enter position [x y]: ');
    ii=round(r(1)/h + 1/2) % Place charge at nearest grid point
    jj=round(r(2)/h + 1/2);
    q = input('Enter charge density: ');
    rho(ii,jj) = rho(ii,jj)+ q/h^2;
end


%compute P (matrix)

cx = cos((2*pi/N)*(0:N-1));
cy = cx;
numerator = -h^2/(2*eps0);
tinynumber = 10^(-20); %avoid division by 0
for i=1:N
    for j=1:N
        P(i,j) = numerator/(cx(i)+cy(j)-2+tinynumber); %alpha = 2pi*sqrt(-1/N)
    end
end

%Compute potential using MFT
rhoT = fft2(rho); %Transfer rho into wavenumber domain
phiT = rhoT.*P; %computing phi in the wavenumber domain
phi = ifft2(phiT); %inv transform phi into coordinates domain;
phi = real(phi); %clean up part due to round-off

%Compute electric field E = -grad phi
[Ex Ey] = gradient(flipud(rot90(phi)));
magnitude = sqrt(Ex.^2 + Ey.^2);
Ex = Ex./magnitude; %normalize
Ey = Ey./magnitude;

%Plot potential and electric field

figure(1); clf;
contour3(x,y,flipud(rot90(phi,1)),35);
xlabel('x'); ylabel('y'); zlabel('\Phi(x,y)');
figure(2); clf;
quiver(x,y,Ex,Ey) %Plot E field with vectors
title('E field (Direction)'); xlabel('x'); ylabel('y');
axis('square'); axis([0 L 0 L]);


