
function [x, u_bar] = solve_newton(l, T, alpha, dx, dt)

% Finite differene (Newton's optimization) method for GDIM thin film 
% equation with periodic boundary condition.
%
% Inputs
%   l: Simulation domain radius, eg. \Omega=[-l,l]. Default value is 2.
%   T: Final time of the simulation. Default value is 50. 
%   alpha: Kernel width. The small scale parameters of the GDIM-TFE. 
%       Default value is 0.05.
%   dx: Grid spacing. Default value is 0.02.
%   dt: Step size of time marching. Default value is 0.01 (5000 steps). 
%
% Outputs 
%   x: 1D array. Spacial grid of the simulation.
%   u_bar: 1D array. Solution of PDE at final time. 
%
% Dependencies
%   my_diag.m

arguments 
    l double = 2;
    T double = 50;
    alpha double = 0.05;
    dx double = 0.02;
    dt double = 0.01;
end

x = -l:dx:l-dx;
Nx = length(x);
t = 0:dt:T;
Nt = length(t);
iter = 20;          % nls maximum interations
tol = 1e-9;         % nls break tolerance

% operators
I = eye(Nx);
D1 = my_diag([-1,0,1], Nx, 'periodic') / (2*dx);
D2 = my_diag([1,-2,1], Nx, 'periodic') / dx^2;
D3 = my_diag([-1,2,0,-2,1], Nx, 'periodic') / (2*dx^3);
kernel = 1/(4*alpha^2) * (alpha + abs(x)) .* exp(-abs(x)/alpha);
K = my_diag(kernel, Nx, 'periodic') * dx;
K = full(K);
Kinv = (I-alpha^2*D2)^2;
dtD1_bar = dt*K*D1;

% initial condition
u = 1.5 * (0.5^2 - x'.^2);
u(abs(x')>=0.5) = 0;
u_bar = K * u;

for i = 2:Nt
    
    v_bar = u_bar;
    
    for j = 1:iter
        
        v = Kinv * v_bar;
        mob = v .* v_bar.^2;
        dkappa = D3 * v_bar;
        Jmob = v_bar.^2 .* Kinv + (2*v_bar.*v) .* I;
        
        % objective function and jacobian
        F = v_bar + dtD1_bar * (mob.*dkappa) - u_bar;
        J = I + dtD1_bar * (Jmob.*dkappa + mob.*D3);
        
        % descent direction
        dv_bar = -J \ F;
        [a, res] = fminbnd(@(a) cost(v_bar+a*dv_bar, u_bar), 0, 1);
        v_bar = v_bar + a*dv_bar;
        if res < tol
            break
        end

    end

    u_bar = v_bar;
    
    if mod(i,10) == 0
        plot(x, u_bar, x, v);
        title(sprintf('t = %.2f, iter = %d, res = %.2e', t(i), j, res))
        drawnow
    end
        
end

du_bar = D1 * u_bar;

function val = cost(v_bar, u_bar)
    
    f = v_bar + dtD1_bar*((Kinv * v_bar) .* v_bar.^2 .* (D3*v_bar)) - u_bar;
    val = norm(f);

end

end