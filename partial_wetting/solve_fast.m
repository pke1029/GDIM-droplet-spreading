
function [m, v, u, du] = solve_fast(l, dx, T, alpha, A0, r0, chi)

% Particle method with fast summation algorithm. (Partial wetting)
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
%       ode solver. Default value is 10.
%   alpha: Kernel width. The small scale parameters of the GDIM-TFE. 
%       Default value is 0.05.
%   dx: Grid spacing. For the particle method, this determine the initial
%       spacing between particles. Default value is 0.02 (200 particles).
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
Nx = length(x);

% initial condition
a = clamp(x-dx/2, -r0, r0);
b = clamp(x+dx/2, -r0, r0);
m = A0*3/(4*r0)*((b-a)-1/(3*r0^2)*(b.^3-a.^3))';

% kernel
N = 1/(4*alpha^2);
kernel = @(x) N * (alpha + abs(x)) .* exp(-abs(x)/alpha);
dkernel = @(x) -1/(4*alpha^3) * x .* exp(-abs(x)/alpha);

% solve ODEs
options = odeset('InitialStep', 0.0005);
[t, X] = ode15s(@(t,x) ODE(m, x), [0, T], x', options);

v = X(end, :);
m = m';
u = sum(m .* kernel(v-v'), 2);
du = sum(m .* dkernel(v-v'), 2);

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
    ell = sqrt(sum(w .* h_bar));
    dh_bar = -N/alpha*(em.*(x.*a1-a2) + ep.*(x.*b1-b2));    

    dPi = -2*chi*A0^2/ell^4 * dh_bar;
    dxdt = h_bar.^2 .* (dkappa - dPi);

end

end