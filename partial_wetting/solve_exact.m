
function [x, h, h_bar, A, B, C, D, xi_star, r_star] = solve_exact(l, dx, alpha, A0, chi)

% Solve for the exact solution of the equilibrium droplet profile.
%
%   Solve boundary value problem for the constant of integration A, B, C, 
%   and D of the equilibrium profile. For discription, see Section 6. 
%
% Inputs
%   l: Domain of evaluation. Default value is 2.
%   dx: Grid spacing of evaluation. Default value is 0.02.
%   alpha: Kernel width. The small scale parameters of the GDIM-TFE. 
%       Default value is 0.05.
%   A0: Area/mass of the droplet. Default value is 1.
%   chi: Droplet width correction factor. chi should be close to 1. Default 
%       value is 1.1602. 
%
% Output
%   x: 1D array. Spacial grid of the constructed solution.
%   h: 1D array. Noisy surface height of the droplet profile.
%   h_bar: 1D array. Smoothen surface height of the droplet profile.
%   A, B, C, D: double. Coefficient of the equilibrium solution. 
%   xi_star, r_star: double. Wave number and droplet core radius of the
%       equilibrium solution.
%
% Dependencies
%   -

arguments 
    l double = 2;
    dx double = 0.02;
    alpha double = 0.05;
    A0 double = 1;
    chi double = 1.1602;
end

% define symbolic variables
xi_sym = sym('xi_sym');
r_sym = sym('r_sym');
y = sym('y');

% linear system
M = [1/xi_sym*sin(xi_sym*r_sym), r_sym, 0, 0;
     1/(1+alpha^2*xi_sym^2)^2/xi_sym*sin(xi_sym*r_sym), r_sym, alpha*exp(-r_sym/alpha), alpha*(r_sym+alpha)*exp(-r_sym/alpha);
     1/(1+alpha^2*xi_sym^2)^2*cos(xi_sym*r_sym), 1, -exp(-r_sym/alpha), -r_sym*exp(-r_sym/alpha);
     -1/(1+alpha^2*xi_sym^2)^2*xi_sym*sin(xi_sym*r_sym), 0, 1/alpha*exp(-r_sym/alpha), + (r_sym/alpha-1)*exp(-r_sym/alpha)];
RHS = [A0/2;
       A0/2;
       0;
       0];
% solve coefficient 
coef = M \ RHS;
A = coef(1);
B = coef(2);
C = coef(3);
D = coef(4);

% inner product
h = A*cos(xi_sym*y) + B;
h_bar = A/(1+alpha^2*xi_sym^2)^2 * cos(xi_sym*y) + B;
ell_sym = 2*int(h*h_bar, y, 0, r_sym);
ell = @(xi,r) double(subs(ell_sym, {xi_sym, r_sym}, {xi, r}));

% equations
eq1_sym = -A/(1+alpha^2*xi_sym^2)^2*xi_sym^2*cos(xi_sym*r_sym) - (C/alpha^2-2*D/alpha+D/alpha^2*r_sym)*exp(-r_sym/alpha);
eq1 = @(xi,r) double(subs(eq1_sym, {xi_sym, r_sym}, {xi, r}));
eq2 = @(xi,r) xi^2 - 2*chi*A0^2/ell(xi,r)^2;

% solution
xi_star = solve_eq2_along_curve();
r_star = solve_eq1(xi_star);

% substitude value
A = double(subs(A, {xi_sym, r_sym}, {xi_star, r_star}));
B = double(subs(B, {xi_sym, r_sym}, {xi_star, r_star}));
C = double(subs(C, {xi_sym, r_sym}, {xi_star, r_star}));
D = double(subs(D, {xi_sym, r_sym}, {xi_star, r_star}));

% construct solution
x = my_centered_array(l, dx);
h = A*cos(xi_star * x') + B;
h(abs(x)>r_star) = 0;

% construct h_bar
h_bar = zeros(length(x), 1);
core = abs(x)<r_star;
h_bar(core) = A/(1+alpha^2*xi_star^2)^2 * cos(xi_star*x(core)) + B;
h_bar(~core) = C*exp(-abs(x(~core))/alpha) + D*abs(x(~core)).*exp(-abs(x(~core))/alpha);

function r = solve_eq1(xi)
    
    r_guess = pi/xi;
    r = fzero(@(r) eq1(xi, r), r_guess);

end

function val = get_eq2_val(xi)

    r = solve_eq1(xi);
    val = eq2(xi, r);

end

function [xi_star] = solve_eq2_along_curve()

    xi_guess = 3;
    xi_star = fzero(@(xi) get_eq2_val(xi), xi_guess);

end

end