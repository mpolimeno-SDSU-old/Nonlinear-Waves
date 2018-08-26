function psi = poisson_attempt3(f,N)
L = 1;
Dd1 = 0:N/2-1;
Dd1(1,1) = 10^(-6);
Dd2 = -N/2:-1;
Dds = 2*1i*pi/2*[Dd1 Dd2]';
Dy = kron(Dds,ones(N,1));
Dx = kron(ones(N,1),Dds);

Dy2 = Dy.^2;
Dx2 = Dx.^2;

Lop = Dy2+Dx2;

% [s1,s2] = size(Lop);
% eps = 1e-10;

% invLop =-1./(Lop);

f = fft2(f);
f = reshape(f',N^2,1);

psi = f./Lop;

psi = reshape(psi.',N,N).';
% psi(1,1) = 0;
psi = real(ifft2(psi));

end