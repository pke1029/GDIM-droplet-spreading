
function A = my_diag(x, N, mode)

% Generate square banded diagonal matrix

if nargin < 3, mode = 'clip'; end

p = length(x);

if strcmp(mode, 'periodic')
    
    % periodic banded matrix
    A = zeros(N);
    y = zeros(1, N);
    y(1:p) = x;
    y = circshift(y, -floor(p/2));
    for i = 1:N
        A(i,:) = circshift(y, i-1);
    end
    
else
    
    % banded diagonal
    A = spdiags(repmat(x, N, 1), (1:p)-floor(p/2)-1, N, N);
    
end

end