
function [v, m, u, du] = solve_particle(l, T, alpha, dx)

% Particle method for GDIM thin film equation with periodic boundary
% condition.
%
% Inputs
%   l: Simulation domain radius, eg. \Omega=[-l,l]. For the particle 
%       method, this is the domain which the particles are initialised, 
%       but the particles are allowed to move outside this 
%       domain during evolution. Default value is 2.
%   T: Final time of the simulation. The time marching is decided by the 
%       ode solver. Default value is 50.
%   alpha: Kernel width. The small scale parameters of the GDIM-TFE. 
%       Default value is 0.05.
%   dx: Grid spacing. For the particle method, this determine the initial
%       spacing between particles. Default value is 0.02 (200 particles).
% 
% Outputs 
%   v: 1D array. Particle positions at final time. 
%   m: 1D array. Particle weights. 
%   u,du: 1D array. Droplet surface height and slope at final time. 
%       These should be plotted against v. eg. plot(v, u, '-o')
%
% Dependencies
%   my_centered_array.m

arguments 
    l double = 2;
    T double = 50;
    alpha double = 0.05;
    dx double = 0.02;
end

x = my_centered_array(l, dx);
u0 = 1.5 * (0.5^2 - x.^2);
u0(abs(x')>=0.5) = 0;
m = u0 * dx;

% kernels
kernel = @(x) 1/(4*alpha^2) * (alpha + abs(x)) .* exp(-abs(x)/alpha);
dkernel = @(x) -1/(4*alpha^3) * x .* exp(-abs(x)/alpha);
dddkernel = @(x) 1/(4*alpha^4) * (2*sign(x) - x/alpha) .* exp(-abs(x)/alpha);

% solve ode
[~, V] = ode45(@(t,v) ODE(m, v), [0, T], x');

v = V(end,:);
xx = v-v';
xx = mod(xx+l, 2*l)-l;
u = sum(m' .* kernel(xx), 1);
du = sum(m' .* dkernel(xx), 1);

function dvdt = ODE(w, v)

    vv = v-v';
    vv = mod(vv+l, 2*l)-l;                  % periodic boundary condition
    
    h_bar = sum(w .* kernel(vv), 2);        % regularize solution
    mu = h_bar.^2;                          % mobility
    dkappa = sum(w .* dddkernel(vv), 2);    % gradient of mean curvature
    
    dvdt = mu .* dkappa;
    
end

end