
% Particle method for Geometric Diffuse Interface Method (GDIM) Thin film 
% equation with infinite boundary condition
%
%   This is the same implementaion as solve_particle.m but with code for
%   generating plots seen in the paper. 
%
% Parameters
%   l: Simulation domain radius, eg. \Omega=[-l,l]. For the particle 
%       method, this is the domain which the particles are initialised, 
%       but the particles are allowed to move outside this 
%       domain during evolution. Default value is 2.
%   T: Final time of the simulation. The time marching is decided by the 
%       ode solver. Default value is 50.
%   alpha: Kernel width. The small scale parameters of the GDIM-TFE. 
%       Default value is 0.05.
%   dx: Grid spacing. For the particle method, this determine the initial
%       spacing between particles. Default value is 0.02 (200 particles),
%       use 0.005 for more accurate results. 
%
% Dependencies
%   clamp.m

% parameters
dx = 0.02;                  % grid spacing
l = 2;                      % domain size
N = 2*l/dx + 1;             % resolution 
x = linspace(-l, l, N);     % space grid
T = 50;                     % final time
alpha = 0.05;               % length scale

% initial condition
A = 1/4;                                        % area
r = 0.5;                                        % radius
a = clamp(x-dx/2, -r, r);
b = clamp(x+dx/2, -r, r);
m = A*3/(4*r)*((b-a)-1/(3*r^2)*(b.^3-a.^3));   % particle weights

% solve ODEs
opts = odeset('Stats', 'on');
[t, V] = ode15s(@(t,v) ODE(m, alpha, v), [0, T], x', opts);

% --------------------------- plots ---------------------------------------

% 1. space time plot
U = sum(permute(m, [1,3,2]) .* kernel(V - permute(V, [1,3,2]), alpha), 3);
subplot(2,3,1);
[M, T] = meshgrid(m, t);
mesh(T, V, U);
xlim([0, 50]);

% 2. solution at time T
v = V(end,:);
u = sum(m' .* kernel(v-v', alpha), 1);
du = sum(m' .* dkernel(v-v', alpha), 1);
subplot(2,3,2);
plot(v, u, '-o', v, du, '-o');

% 3. similarity variable
subplot(2,3,3);
contourf(V.*T.^(-1/7), T, U.*T.^(1/7))
xlim([-1, 1])

% 4. similarity solution
f0 = [t(end)^(1/7)*U(end,floor((N+1)/2)); 0; -1.15];
[eta, f] = ode45(@(t,f) [f(2); f(3); t/7 * f(1)^(-2)], [0, 1.5], f0);
subplot(2,3,4);
plot(t(end)^(-1/7)*V(end,:), t(end)^(1/7)*U(end,:), eta, f(:,1), '-o')
ylim([0, 0.3])

% 5. trajectory/characteristic
subplot(2,3,5);
mesh(T, V, M)
view([0, 90]);
xlim([0, 50]);
ylim([0,1.3]);

% 6. tanner's law
subplot(2,3,6);
mesh(T, V, M)
set(gca, 'view', [0, 90], 'xlim', [t(2), t(end)], 'ylim', [0.4, 1.2], 'XScale', 'log', 'YScale', 'log')
hold on
plt = plot3(t, t.^(1/7)/1.51, t*0+1, '-r', 'linewidth', 2, 'displayname', 't^{1/7}');
hold off
legend([plt])

% --------------------------- functions -----------------------------------

function val = kernel(x, alpha)
val = 1/(4*alpha^2) * (alpha + abs(x)) .* exp(-abs(x)/alpha);
end

function val = dkernel(x, alpha)
val = -1/(4*alpha^3) * x .* exp(-abs(x)/alpha);
end

function val = dddkernel(x, alpha)
val = 1/(4*alpha^4) * (2*sign(x) - x/alpha) .* exp(-abs(x)/alpha);
end

function dvdt = ODE(w, alpha, v)

V = v-v';

h_bar = sum(w .* kernel(V, alpha), 2);          % regularize solution
mu = h_bar.^2;                                  % mobility
dkappa = sum(w .* dddkernel(V, alpha), 2);      % gradient of mean curvature

dvdt = mu .* dkappa;

end
