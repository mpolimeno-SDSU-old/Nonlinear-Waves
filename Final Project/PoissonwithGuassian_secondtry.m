% Nhat Hoang April 17th 
% Solve the Poisson Equation using SOR method
% Solve the 

clear; close all;
% part 1L initalize variabes:
L = 1;
dx = .025;
dy = dx;
h= dx;
n = round(2*L/dx);
m = n;
x = linspace(-L,L,n)';
y = linspace(-L,L,m)';
T_0 = zeros(n,m);

bb = 2/(1 + pi*h); % this is my 

% Boundary Condition
% At x= 0; 
for i =1:n
    T_0(i,1) = 0;  
% At x = 1
    T_0(i,m) = 0;
end
% At y = 0;
for j =1:m
    T_0(1,j) = 0;
% At y = 1;
    T_0(m,j) = 0;
end


T= T_0;

% Define source terms
w_0 = zeros(n,m);  % Original vorticity

for i = 1:n
   for j=1:m
         w_0(i,j)= exp(-(x(i)^2*2 + y(j)^2/20));
    end
end


figure(1);
surf(x,y,T_0);
title('Initial and Boundary Condition');
xlabel('x');
ylabel('y');
zlabel('z');

% plot actual solution to compare later


max_step = 400;
%Lnorm = 100;
%l = 0;


for l=1:max_step
 for i=2:n-1
     % Update new vorticity:
     w= w_0;      
     % Boundary Conditions
     T_old = T;
     for j=2:m-1
        T(i,j)=bb*0.25*(T(i+1,j)+...
        T(i,j+1)+T(i-1,j)+T(i,j-1)-h^2*w(i,j))+(1.0-bb)*T(i,j);
        % Boundary Conditions
       
     end
 end
 % find residual
    res = (T-T_old);
    res = reshape(res,1,[]); % flatten matrix, [] means maintain the same number of elements
    res = norm(res);
    if res < 2*10^-6
        break
    end
 
 if rem(l,10) == 0
     figure(3);
     surf(x,y,T);
     title([' SOR for Poisson Equation, at iteration ',num2str(l),' with'...
        ' h =',num2str(h)]);
     xlabel('x')
     ylabel('y')
     zlabel('z')
     %axis([ -1 1 -1 1 -10 10]);
     getframe(gcf);
     pause(.1);
 end
 
 
end

% Update 



figure(4);
contour(T)
title(['Contour Plot, residual value ', num2str(res), 'with h = ',num2str(h)]);
