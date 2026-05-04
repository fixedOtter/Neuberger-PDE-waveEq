% newtons method for PDE's
% written by fixedotter

clear; clc;

% boring const
n = 20; dx = 1/n; n_half = round(n/2);
ds = 0.1; % s_f = 3*pi^2; % final s
s_f = 0;

% array of s 
s_array = 0:ds:s_f;
u_array = zeros(1,length(s_array)); % to store solutions for each s

% basically just the grid
x = linspace(0,1,n+1)'; xs = x(2:end-1); % full and interior pts in interval 
[YS,XS,ZS] = meshgrid(xs,xs,xs); % grid of interior points

% defining the laplacian in 2d
I  = speye(n-1);
I_sq = speye((n-1)^2);
D2 = spdiags(kron(ones(n-1,1),[1 -2 1]),-1:1,n-1,n-1)/dx^2;

% Laplacian = D_2x + D_2^y
L  = kron(I,D2) + kron(D2, I);

% laplacian w/ cube
L_cub = kron(I_sq,D2) + kron(L,I);

% compute first 6 eigen pairs
[V,D] = eigs(-L_cub,6,'sm'); 
disp(D)
% 6 smallest DIFFERENT eigenvalues are 3*pi^2, 6*pi^2, 9*pi^2, 11*pi^2, 12*pi^2, 14*pi^2
% since there's a linear combination, we basically get 3 of each of them after 3*pi^2, which is what we see in the diag(D) output

% function & it's jacobian with s to move along solution space?
F = @(s,u) L_cub*u + s*u + u.^3; % actual funct with s to move along the solutions
J = @(s,u) L_cub + spdiags(s + 3*u.^2, 0, (n-1)^3, (n-1)^3); % jacobian with s

% for some reason 8 times the whole thing is a good guess but 9 or 7 isnt
f = @(x,y,z) sin(pi*x).*sin(pi*y).*sin(pi*z)*8; % init guess


% putting her into a column for easy matrix ops
u = reshape(f(XS,YS,ZS), (n-1)^3, 1);

% % do again with s & see if we can loop over the solutions to plot a grid
% for s = 0:ds:2*pi^2
%   chi = J(s,u)\F(s,u); % next newton step
%   u = u - chi; % update
%   if (norm(chi) < 1e-6) % stopping condition
%     break
%   end
% end

% just doing s = 0 for now
for s = 0:ds:s_f
  % newtons meth step
  for i = 1:250
    chi = J(s,u)\F(s,u); % next newton step
    u = u - chi; % update
    fprintf("iteration %d, norm of chi: %e\n", i, norm(chi));
    disp(max(abs(u))); % max value of u
    if (norm(chi) < 1e-6) % stopping condition
      break
    end
  end
  u_array(s+1) = max(abs(u)); % store the max value of u for this s
end

% bringing it back to 3d for plotting
u = reshape(u, n-1, n-1, n-1); % reshape function
Sol_3d = zeros(n+1,n+1,n+1); % boundary conditions
Sol_3d(2:end-1,2:end-1,2:end-1) = u; % feeding in the solution
[X,Y,Z] = meshgrid(x,x,x); % grid for plotting


% can we save this to a mf gif
gifFile = "neubergers_cube_w_newtonsMeth.gif";
% for k = 0:n
%   figure("Visible",'off')
%   colormap(jet); % shading interp
%   Z_slice = Z(:,:,k+1); % slicing the solution at z = k*ds
%   zlim([-1,1]); % keeping the z axis fixed so it's less confusing
%   xlim([0,1]); ylim([0,1]); % keeping x and y axis fixed too
%   surf(X(:,:,k+1),Y(:,:,k+1),Z_slice,'edgecolor','none'); % plotting the slice of the solution
%   titleName = strcat('z = ', num2str(k*dx));
%   title(titleName)
%   view(42,-42);
%   % exportgraphics(gcf, strcat('neubergers_cube_w_newtonsMeth_s=', num2str(k*ds), '.png')); % saving the figure
%   exportgraphics(gcf,gifFile,Append=true); % saving to the gif
% end

% surface figure
figure
surf(X(:,:,n_half),Y(:,:,n_half),Sol_3d(:,:,n_half),'edgecolor','none'); % plotting the slice of the solution
title('Slice of the solution at z = 0.5')
% zlim([0 10E-34])
view(42,-42);

scatX = reshape(X, [],1,1);
scatY = reshape(Y, [],1,1);
scatZ = reshape(Z, [],1,1);
scatSol_3d = reshape(Sol_3d, [],1,1);

% scatter3d figure
figure
scatter3(scatX,scatY,scatZ)
colormap(gca,"parula")

figure
plot(s_array, u_array, '-o');
xlabel('s value');
ylabel('max absolute value u');
title('max |u| vs s');
grid on;

% figure
% colormap(jet); % shading interp
% isosurface(X,Y,Z,0); % plotting the isosurface of the solution
% title('Isosurface of the solution at s = 3*pi^2')
% view(42,-42);



  % what chat said
  % figure('Visible','off','Name','gif');
  % colormap(jet); % shading interp
  % isosurface(X,Y,Z,k*ds); % plotting the isosurface of the solution
  % titleVar = strcat('s = ', num2str(k*ds));
  % title(titleVar)
  % view(42,-42);
  % saveas(gcf, strcat('wave_sol_s=', num2str(k*ds), '.png'));

% isosurface(X,Y,Z,0); % plotting the isosurface of the solution