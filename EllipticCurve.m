clear; clc;

x = 0:70;
status = zeros(size(x));
status_zero = zeros(size(x));
y_square = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    y_square(i) = mod(xi^3+xi+28,71);
    if SquareAndMultiply(y_square(i),35,71) == 1
        status(i) = 1;
    elseif SquareAndMultiply(y_square(i),35,71) == 0
        status_zero(i) = 1;
    end
end
Nb_pts = nnz(status)*2+nnz(status_zero)+1

function z = SquareAndMultiply(x,c,n)
% Computes modular exponentiation x^c mod n
    ci = de2bi(c);
    z = 1;
    for i = length(ci):-1:1
        z = mod(z^2,n);
        if ci(i) == 1
            z = mod(z*x,n);
        end
    end
end