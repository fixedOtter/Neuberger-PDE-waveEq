% written by gunnar march 19 2026

% what are we solving? the wave equation u_tt = delta*u, with 0-D BC and initial conditions

% first we do the laplacian
n = 1024; dx = 1/n;
x = linspace(0,1,n+1)'; xs = x(2:end-1); % full and interior pts in interval 
[YS,XS] = meshgrid(xs,xs); % grid of interior points

% we need to get the L matix 
I  = speye(n-1);
D2 = spdiags(kron(ones(n-1,1),[1 -2 1]),-1:1,n-1,n-1)/dx^2;
L  = kron(I,D2) + kron(D2, I);   % Laplacian = D_2^x + D_2^y

% RHS of the ODE system
F = @(t,w) [w((n-1)^2+1:end); L * w(1:(n-1)^2)];

%  initial conditions
% position
% pure mode
% f = @(x,y) sin(pi*x).*sin(pi*y); % pos
% f = @(x,y) sin(3*pi*x).*sin(3*pi*y); % pos
f = @(x,y) (x>.4).*(x<.6).*(y>.4).*(y<.6); % pos
% velocity
g = @(x,y) zeros(size(x)); % vel

% initial conditons into w0
% really F(xs,ys) but ys = xs since we are on a square
w0 = [reshape(f(XS,YS), (n-1)^2,1); reshape(g(XS,YS), (n-1)^2,1)]; % initial conditions

% time span
tf = 4;
m = 1;
t = linspace(0,tf,m+1);

% solve the ODE system
[t,w] = ode45(F, t, w0);

% reminder: w = [u; v] so we need to reshape the first half of w to get u

% plotting the solution
% i wanna make a movie
gifFile = "wave.gif";
z = zeros(n+1);
for i=1:length(t)
  % reshaping the solution so we get u
  z(2:end-1,2:end-1) = reshape(w(i,1:(n-1)^2), n-1, n-1);
  % plotting the fig
  figure('Visible','on','Name','gif');
  colormap(jet); % shading interp;
  surf(x,x,z,'edgecolor','none'); 
  % keeping the z axis fixed so it's less confusing
  zlim([-1,1]);
  titleVar = strcat('t = ', num2str(t(i)));
  % title(sym(strcat('Solving U_tt =', num2str(t(i))))); 
  title(titleVar)
  clim([-1,1]);
  view(42,-42);

  % saving to the video
  % frame = getframe(gcf);
  exportgraphics(gcf,gifFile,Append=true);
end