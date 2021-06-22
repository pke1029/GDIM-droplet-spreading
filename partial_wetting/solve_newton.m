
function [x, v, u_bar, D1, D2, D3] = solve_newton(l, dx, T, alpha, A0, r0, chi, dt)

% Finite differene (Newton's optimization) method for GDIM thin film 
% equation with periodic boundary condition. (Partial wetting)
%
% Inputs
%   l: Simulation domain radius, eg. \Omega=[-l,l]. Default value is 2.
%   T: Final time of the simulation. Default value is 10. 
%   alpha: Kernel width. The small scale parameters of the GDIM-TFE. 
%       Default value is 0.05.
%   dx: Grid spacing. Default value is 0.02.
%   A0: Area/mass of the initial droplet profile. Default value is 1.
%   r0: Radius of the initial droplet base. Default value is 0.5. 
%   chi: Droplet width correction factor. chi should be close to 1. Default 
%       value is 1.1602. 
%   dt: Step size of time marching. Default value is 0.01 (5000 steps). 
%
% Outputs 
%   x: 1D array. Spacial grid of the simulation.
%   u_bar: 1D array. Solution of PDE at final time. 
%
% Dependencies
%   my-centered_array.m
%   my_diag.m

arguments 
    l double = 2;
    dx double = 0.02;
    T double = 10;
    alpha double = 0.05;
    A0 double = 1;
    r0 double = 0.5;
    chi double = 1.1602;
    dt double = 0.01;
end

x = my_centered_array(l, dx);
Nx = length(x);
t = 0:dt:T;
Nt = length(t);
iter = 20;          % nls maximum interations
tol = 1e-7;         % nls break tolerance

% operators
I = eye(Nx);
D1 = my_diag([-1,0,1], Nx, 'periodic') / (2*dx);
D2 = my_diag([1,-2,1], Nx, 'periodic') / dx^2;
D3 = my_diag([-1,2,0,-2,1], Nx, 'periodic') / (2*dx^3);
kernel = 1/(4*alpha^2) * (alpha + abs(x)) .* exp(-abs(x)/alpha);
K = my_diag(kernel, Nx, 'periodic') * dx;
K = full(K);
Kinv = (I-alpha^2*D2)^2;
D1_bar = K*D1;
dtD1_bar = dt*D1_bar;

% initial condition
u = A0*3/(4*r0)*(1-(x'/r0).^2);
u(u<0) = 0;
u_bar = K * u;

% potential function
Pi = @(h, h_bar) -2*chi*A0^2/(sum(h.*h_bar)*dx)^2 .* h_bar;

% mobility function
mob = @(h, h_bar) h.*h_bar.^2;
grad_mob = @(h, h_bar) h_bar.^2.*Kinv + 2*h.*h_bar.*I;

for i = 2:Nt
    
    u = Kinv * u_bar;
    v_bar = u_bar;
    
    for j = 1:iter
        
        v = Kinv * v_bar;
        mu = mob(v, v_bar);
        dkappa = D1 * (D2*v_bar - Pi(u, u_bar));
        
        % decent direction
        F = v_bar + dtD1_bar * (mu.*dkappa) - u_bar;
        J = I + dtD1_bar * (grad_mob(v, v_bar).*dkappa + mu.*D3);
        dv_bar = -J \ F;
        
        % find optimum stepsize
        [a, res] = fminbnd(@(a) cost(v_bar+a*dv_bar, u, u_bar), 0, 1);
        v_bar = v_bar + a*dv_bar;
        if res < tol
            break
        end

    end
    
    u_bar = v_bar;
    
    if mod(i, 100) == 2
        plot(x, u_bar, x, v);
        title(sprintf('t = %.2f', t(i)))
        drawnow
    end
        
end

du_bar = D1 * u_bar;

function val = cost(v_bar, u, u_bar)
    
    f = v_bar + dtD1_bar*(mob(Kinv*v_bar, v_bar) .* (D1*(D2*v_bar - Pi(u, u_bar)))) - u_bar;
    val = norm(f);

end

end
