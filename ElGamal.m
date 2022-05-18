clear; clc;

p = [1 2 0 1];  % p = x^3+2x^2+1
c = 11;  % exponent

% Examples of y1 and y2
y1 = [2 2];  % y1 = 2x+2 (letter H)
y2 = [2 0 0];   % y2 = 2x^2 (letter R)

% answer = PolyReduction(PolyMultiplication(y2,PolySquareAndMultiply(PolyInverse(y1,p),c,p)),p)
answer = PolyReduction(PolyMultiplication(y2,PolyInverse(PolySquareAndMultiply(y1,c,p),p)),p)

function z = PolySquareAndMultiply(f,c,p)
% Square and multiply algorithm
    ci = de2bi(c);
    z = 1;
    for i = length(ci):-1:1
        z = PolyReduction(PolyMultiplication(z,z),p);
        if ci(i) == 1
            z = PolyReduction(PolyMultiplication(z,f),p);
        end
    end
end

function f_inv = PolyInverse(f,p)
    for i = 0:2
        for j = 0:2
            for k = 0:2
                if not(i==0&&j==0&&k==0)
                    candidate = [i j k];
                    if length(PolyReduction(PolyMultiplication(f,candidate),p))==1 && PolyReduction(PolyMultiplication(f,candidate),p)==1
                        f_inv = candidate;
                    end
                end
            end
        end
    end
end

function result = PolyReduction(f,p)
% result = f(x) mod p(x) in Z_3[x]
    while length(f)>=length(p)
        f = [EquDegPolysubtraction(f(1:length(p)),p),f(length(p)+1:end)];
    end
    result = f;
end

function result = EquDegPolysubtraction(f,g)
% Equal degree polynomial subtraction in Z_3[x]
    if f(1)~=g(1)
        result = mod(f-2*g,3);
    else
        result = mod(f-g,3);
    end
    result = unpad(result);
end

function vec = unpad(w)
%   Takes a row vector, removes all zeros at the beginning of the vector,
%   and returns the remaining elements of the vector
    index = find(w ~= 0, 1, 'first');
    vec = w(index:end);
end

function result = PolyMultiplication(f,g)
% f(x)*g(x) in Z_3[x]
    result = mod(conv(f,g),3);
end