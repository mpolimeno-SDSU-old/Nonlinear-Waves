function psi = poisson_attempt2(f,N)

% f = fft2(f);
f = f(:);
% f = [0; f(:); 0];
[D,x] = cheb(N-1);
D(N,:) = zeros(1,N);
D(1,:) = zeros(1,N);
% D(:,N) = zeros(1,N);
% D(:,1) = zeros(1,N);
D2 = D^2; 
y = x;

% D2 = D2(1:N,1:N); I = eye(N);

I = eye(length(D2));

L = kron(I,D2)+kron(D2,I);

% L = L(2:N,2:N);

% size(L)
% size(f)


psi = f\L;

psi = reshape(psi,N,N);
psi(N,:) = zeros(1,N);
psi(1,:) = zeros(1,N);
psi(:,N) = zeros(1,N);
psi(:,1) = zeros(1,N);

% psi = real(ifft2(psi));


end