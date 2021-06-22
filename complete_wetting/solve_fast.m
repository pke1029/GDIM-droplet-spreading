
function [v, m, u, du] = solve_fast(l, T, alpha, dx)

% Particle method with fast summation algorithm
%
%   The fast summation algorithm reduces the computational complexity of
%   the particle method from O(N^2) to O(N) where N is the number of
%   particles. This is achieved by using a recursive relation to evaluate
%   (components of) h_bar and dkappa. For the derivation, please see
%   Appendix A. 
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
Nx = length(x);
u0 = 1.5 * (0.5^2 - x'.^2);
u0(abs(x')>=0.5) = 0;
m = u0 * dx;

% kernels
N = 1/(4*alpha^2);
kernel = @(x) 1/(4*alpha^2) * (alpha + abs(x)) .* exp(-abs(x)/alpha);
dkernel = @(x) -1/(4*alpha^3) * x .* exp(-abs(x)/alpha);
dddkernel = @(x) 1/(4*alpha^4) * (2*sign(x) - x/alpha) .* exp(-abs(x)/alpha);

% solve ode
options = odeset('InitialStep', 0.001);
[~, V] = ode45(@(t,v) ODE(m, v), [0, T], x', options);

v = V(end,:);
m = m';
u = sum(m' .* kernel(v-v'), 1);
du = sum(m' .* dkernel(v-v'), 1);

function dxdt = ODE(w, x)
    
    ep = exp(x/alpha);
    em = exp(-x/alpha);
    
    a1 = zeros(Nx, 1);
    a2 = zeros(Nx, 1);
    b1 = zeros(Nx, 1);
    b2 = zeros(Nx, 1);
    
    for i = 1:Nx-1
    
        a1(i+1) = a1(i) + w(i)*ep(i);
        a2(i+1) = a2(i) + w(i)*x(i)*ep(i);

        j = Nx-i;
        b1(j) = b1(j+1) + w(j+1)*em(j+1);
        b2(j) = b2(j+1) + w(j+1)*x(j+1)*em(j+1);

    end
    
    h_bar = N*em.*((alpha+x).*a1-a2) + w*kernel(0) + N*ep.*((alpha-x).*b1+b2);
    dkappa = N/alpha^2 * (em.*(2-x/alpha).*a1 + 1/alpha*em.*a2 + ep.*(-2-x/alpha).*b1 + 1/alpha*ep.*b2);
    
    dxdt = h_bar.^2 .* dkappa;
    
    % O(N^2) implementation
%     h_bar2 = sum(w' .* kernel(x-x'), 2);
%     dkappa2 = sum(w' .* dddkernel(x-x'), 2);
%     dxdt2 = h_bar2.^2 .* dkappa2;
    
%     plot(x, h_bar)
%     drawnow
    
end

end