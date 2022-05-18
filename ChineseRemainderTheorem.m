clear; clc;

m = [7 11 13];
remainders = [5 3 10];

r = length(m);
M = prod(m);
Mi = M./m;
y = MultiplicativeInverse(Mi,m);
x = mod(sum(remainders.*Mi.*y),M);

function a_inv = MultiplicativeInverse(a,b)
% Computes a_inv mod b using Extended Euclidean Algorithm
% s*a + t*b = r = gcd(a,b)
% Vector inputs possible
    
    n = length(a);
    a_inv = zeros(size(a));
    for i = 1:n
        a0 = a(i); b0 = b(i); t0 = 0; t = 1; s0 = 1;
        s = 0; q = floor(a0/b0); r = a0 - q*b0;
        while r>0
            temp = t0 - q*t; 
            t0 = t;
            t = temp;
            temp = s0 - q*s;
            s0 = s;
            s = temp;
            a0 = b0;
            b0 = r;
            q = floor(a0/b0);
            r = a0 - q*b0;
        end
        r = b0;

        if r == 1
            a_inv(i) = mod(s,b(i));
        end
    end
end