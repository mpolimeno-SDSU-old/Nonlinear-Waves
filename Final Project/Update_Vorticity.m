% Program to Update Vorticity;
% Euler Method

clear; close all;
% part 1: initalize variabes:
L = 2;
dt = .001;
dx = .025;
%lambda = 1/2; % lambda = dt/dx
%dt = lambda*dx;
dy = dx;
h= dx;
n = round(2*L/dx);
m = n;
x = linspace(-L,L,n)';
y = linspace(-L,L,m)';
T_0 = zeros(n,m);
time_max = 5;
max_step = round(time_max/dt);
v = .0001;
time =0;
bb = 2/(1 + pi*h); % this is my SOR coefficient

% Define intial vortice
% Define source terms
w_0 = zeros(n,m);  % Original vorticity

for i = 1:n
   for j=1:m
         w_0(i,j)= exp(-(x(i).^2 + 16.*y(j).^2));
    end
end

w = w_0; % Intialize Vorticity

%Initialize derivatives matrices
dPsi_dx = zeros(n,m);
dPsi_dy = zeros(n,m);
dw_dx = zeros(n,m);
dw_dy = zeros(n,m);
Lap_w = zeros(n,m);
Psi_w_commute = zeros(n,m);

            
tic
for tt = 1:max_step
    Psi = Poisson_func_SOR( w , L ,dx ,bb ); % Solve Poisson Equation by SOR
         
    for i=2:n-1
        for j=2:n-1
            % Define and update derivatives terms:
            dPsi_dx(i,j) = (Psi(i+1,j)+Psi(i-1,j))./(2*dx);
            dPsi_dy(i,j) = (Psi(i,j+1)+Psi(i,j-1))./(2*dy);
            dw_dx(i,j) = (w(i+1,j)+w(i-1,j))./(2*dx);
            dw_dy(i,j) = (w(i,j+1)+w(i,j-1))./(2*dy);
            Lap_w(i,j) = (w(i+1,j)+ w(i-1,j) + w(i,j+1)-w(i,j-1) - 4*w(i,j))/dx^2; % laplacian of w
            Psi_w_commute = dPsi_dx.*dw_dy - dPsi_dy.*dw_dx;
            
            
            %Update the vector through time
            w(i,j) = -w(i,j) + 2*dt.*v.*Lap_w(i,j) - 2*dt*Psi_w_commute(i,j);
        end
    end
    
    time = time + dt;
    %Plot the evolution of the vortice
    figure(1);
    if rem(tt,2) == 0
        pcolor(x,y,w);
        colorbar;
        shading flat; colormap('jet');
        title(['Evolution of the vorticity at time ', num2str(time)]);
    end
%     figure(2);
%         contour(Psi);
%     title('Evolution of the streamfunction')
%   
    
    M(i) = getframe(gcf);
    
end
toc
video =VideoWriter('Finite_Difference_SOR.avi');
open(video);
writeVideo(video,M);
close(video);

