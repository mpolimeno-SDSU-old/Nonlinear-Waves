function w = w0(N,x,y,xx,yy)

xn = length(x);
yn = length(y);

mpx = x(round(xn/2));
mpy = y(round(yn/2));
% h = 0.025;
yl = 40;
xl = 20*yl; %changes width of gaussian hump

a = 1;


w = a*exp(-((xl*(xx).^2)+(yl*(yy).^2)));

% w = a*exp(-((xl*(xx-mpx+.1).^2)+(yl*(yy-mpy+.1).^2))) + a*exp(-((xl*(xx-mpx-.1).^2)+(yl*(yy-mpy-.1).^2)));

% size(w)
%bin's IC
% L = 1;
% dx = .025;
% dy = dx;
% h= dx;
% n = xn;
% m = n;
% 
% w = zeros(n,m); 
% 
% for i = 1:n
%    for j=1:m
%          w(i,j)= exp(-(x(i)^2/h^2+y(j)^2/h^2));
%     end
% end
% w = exp(-xx.^(2/h^2)+yy.^(2/h^2));
% w = exp((-2*(xx).^2 - (yy).^2)/20)

% mesh(w)
% figure()

% size(w)



% w = reshape(w,N,N);

% w = 10*sin(8*xx.*(yy-1));

end
