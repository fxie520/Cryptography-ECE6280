clear;

% Extended Euclidean Algorithm: s*a + t*b = r = gcd(a,b)
% a_inv: the inverse of a modulo b (if exists)

a = 111111111; b = 270153823;

a0 = a; b0 = b; t0 = 0; t = 1; s0 = 1;
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
    a_inv = mod(s,b);
end