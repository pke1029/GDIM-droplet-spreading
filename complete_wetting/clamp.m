
function val = clamp(x, a, b)

% restrict x in the range [a,b]

val = min(max(x, a), b);

end