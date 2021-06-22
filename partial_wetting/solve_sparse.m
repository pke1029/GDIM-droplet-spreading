
function [m, v, u, du] = solve_sparse(l, dx, T, alpha, A0, r0, chi)

% Particle method with sparse implementation. (Partial wetting)
%
%   The sparse implementation merely drop any particles with weight 0 since
%   they do not contribute to the solution. The rest is the same as
%   solve_particle.m
%
% Inputs
%   l: Simulation domain radius, eg. \Omega=[-l,l]. For the particle 
%       method, this is the domain which the particles are initialised, 
%       but the particles are allowed to move outside this 
%       domain during evolution. Default value is 2.
%   dx: Grid spacing. For the particle method, this determine the initial
%       spacing between particles. Default value is 0.02 (200 particles).
%   T: Final time of the simulation. The time marching is decided by the 
%       ode solver. Default value is 10.
%   alpha: Kernel width. The small scale parameters of the GDIM-TFE. 
%       Default value is 0.05.
%   A0: Area/mass of the initial droplet profile. Default value is 1.
%   r0: Radius of the initial droplet base. Default value is 0.5. 
%   chi: Droplet width correction factor. chi should be close to 1. Default 
%       value is 1.1602. 
% 
% Outputs 
%   v: 1D array. Particle positions at final time. 
%   m: 1D array. Particle weights. 
%   u,du: 1D array. Droplet surface height and slope at final time. 
%       These should be plotted against v. eg. plot(v, u, '-o')
%
% Dependencies
%   my_centered_array.m
%   clamp.m

arguments 
    l double = 2;
    dx double = 0.02;
    T double = 10;
    alpha double = 0.05;
    A0 double = 1;
    r0 double = 0.5;
    chi double = 1.1602;
end

x = my_centered_array(l, dx);

% initial condition
a = clamp(x-dx/2, -r0, r0);
b = clamp(x+dx/2, -r0, r0);
m = A0*3/(4*r0)*((b-a)-1/(3*r0^2)*(b.^3-a.^3));

% sparse
x = x(m > 0);
m = m(m > 0);

% kernel
kernel = @(x) 1/(4*alpha^2) * (alpha + abs(x)) .* exp(-abs(x)/alpha);
dkernel = @(x) -1/(4*alpha^3) * x .* exp(-abs(x)/alpha);
dddkernel = @(x) 1/(4*alpha^4) * (2*sign(x) - x/alpha) .* exp(-abs(x)/alpha);

% solve ODEs
[t, X] = ode15s(@(t,v) ODE(m, v), [0, T], x');

v = X(end, :);
u = sum(m .* kernel(v-v'), 2);
du = sum(m .* dkernel(v-v'), 2);

function dvdt = ODE(w, v)

V = v-v';

h_bar = sum(w .* kernel(V), 2);          % regularize solution
mu = h_bar.^2;                           % mobility
dkappa = sum(w .* dddkernel(V), 2);      % gradient of mean curvature

% potential
dh_bar = sum(w .* dkernel(V), 2);
ell = sqrt(w*h_bar);                    % inner product
dPi = -2*chi*A0^2/ell^4 * dh_bar;

dvdt = mu .* (dkappa - dPi);

end

end