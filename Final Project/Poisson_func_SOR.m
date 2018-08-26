function [ Psi ] = Poisson_func( w , L ,dx ,bb ) 

%POISSON_FUNC

%   Solve the Poisson function , input vorticiy w, and the Length L

dy = dx;
h= dx;
n = round(2*L/dx);
m = n;
x = linspace(-L,L,n)';
y = linspace(-L,L,m)';
Psi_0 = zeros(n,m);

% Boundary Condition
% At x= 0; 
for i =1:n
    Psi_0(i,1) = 0;  
% At x = 1
    Psi_0(i,m) = 0;
end
% At y = 0;
for j =1:m
    Psi_0(1,j) = 0;
% At y = 1;
    Psi_0(m,j) = 0;
end


Psi= Psi_0;

% Define source terms
% w_0 = zeros(n,m);  % Original vorticity
% 
% for i = 1:n
%    for j=1:m
%          w_0(i,j)= exp(-(x(i)^2/h^2+y(j)^2/h^2));
%     end
% end


% figure(1);
% surf(x,y,Psi_0);
% title('Boundary Condition');
% xlabel('x');
% ylabel('y');
% zlabel('z');

% plot actual solution to compare later


max_step = 200;
%Lnorm = 100;
%l = 0;


for l=1:max_step
 for i=2:n-1
     % Boundary Conditions
     Psi_old = Psi;
     for j=2:m-1
        Psi(i,j)=bb*0.25*(Psi(i+1,j)+...
        Psi(i,j+1)+Psi(i-1,j)+Psi(i,j-1)-h^2*w(i,j))+(1.0-bb)*Psi(i,j);
        % Boundary Conditions
       
     end
 end
 % find residual
    res = (Psi-Psi_old);
    res = reshape(res,1,[]); % flatten matrix, [] means maintain the same number of elements
    res = norm(res);
    if res < 2*10^-7
        break
    end
 
%  if rem(l,5) == 0
%      figure(3);
%      surf(x,y,Psi);
%      title([' SOR for Poisson Equation, at iteration ',num2str(l),' with'...
%         ' h =',num2str(h)]);
%      xlabel('x')
%      ylabel('y')
%      zlabel('z')
%      %axis([ -1 1 -1 1 -10 10]);
%      getframe(gcf);
%      pause(.1);
%  end
 

% figure(4);
% contour(Psi);
% title(['Contour Plot, residual value ', num2str(res), 'with h = ',num2str(h)]);
Psi;
end

