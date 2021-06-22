
function x = my_centered_array(L, dx)

% Generate array centered at 0

a = dx:dx:L;
x = [-flip(a), 0, a];

end