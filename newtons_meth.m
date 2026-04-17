

% newtons method for PDE's

% defiining our grid
n = 100; dx = 1/n;
x = linspace(0,1,n+1)'; xs = x(2:end-1); % full and interior pts in interval 
[YS,XS] = meshgrid(xs,xs); % grid of interior points

% pulling in the laplacian from other code
I  = speye(n-1);
D2 = spdiags(kron(ones(n-1,1),[1 -2 1]),-1:1,n-1,n-1)/dx^2;
L  = kron(I,D2) + kron(D2, I);   % Laplacian = D_2^x + D_2^y

% anon functions
% the function :)
F = @(u) L*u + u.^3; % actual func
% jacobian
J = @(u) L + spdiags(3*u.^2,0,(n-1)^2,(n-1)^2) ;
% initial guess
f = @(x,y) sin(pi*x).*sin(pi*y)*4; % pos

[Ys,Xs] = meshgrid(xs,xs);
u = reshape(f(Xs,Ys), (n-1)^2,1); % putting in the init

for i = 1:10
  chi = J(u)\F(u); % next newton step
  u = u - chi; % update
  if (norm(chi) < 1e-6) % stopping condition
    break
  end
end

u = reshape(u, n-1, n-1);
Z = zeros(n+1);
Z(2:end-1,2:end-1) = u;
[X,Y] = meshgrid(x,x);
surf(X,Y,Z,'edgecolor','none');