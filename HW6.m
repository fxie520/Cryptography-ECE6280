clear; clc;

% % Problem 2
% mod(SquareAndMultiply(170,49,7879),101)
% mod((52+75*59)*MultiplicativeInverse(49,101),101)
% e1 = mod(52*MultiplicativeInverse(79,101),101);
% e2 = mod(59*MultiplicativeInverse(79,101),101);
% mod(mod(SquareAndMultiply(170,16,7879)*SquareAndMultiply(4567,57,7879),7879),101)

% Problem 3
p = 127; a = 1; b = 26;  % Parameters of the curve
x_A = 2; y_A = 6; InfinityFlag_A = false;  % Point A
x_initial = x_A; y_initial = y_A; InfinityFlag_initial = InfinityFlag_A;
x_current = x_A; y_current = y_A; InfinityFlag_current = false;
for m = 1:53  % Compute B
    [x_current,y_current,InfinityFlag_current] = PointAddition(x_current,y_current,InfinityFlag_current,x_initial,y_initial,InfinityFlag_initial,a,p);
end
x_B = x_current; y_B = y_current; InfinityFlag_B = false;
disp("[x_B,y_B] = ")
disp([x_B,y_B])

x_initial = x_A; y_initial = y_A; InfinityFlag_initial = InfinityFlag_A;
x_current = x_A; y_current = y_A; InfinityFlag_current = false;
for m = 1:74  % Compute kA
    [x_current,y_current,InfinityFlag_current] = PointAddition(x_current,y_current,InfinityFlag_current,x_initial,y_initial,InfinityFlag_initial,a,p);
end
disp("[x_kA,y_kA] = ")
disp([x_current,y_current])

s = mod(MultiplicativeInverse(75,131)*(10+54*88),131)

w = MultiplicativeInverse(s,131)
i = mod(107*10,131)
j = mod(107*88,131)

% Compute iA+jB
x_current = x_A; y_current = y_A; InfinityFlag_current = false;
for m = 1:(i-1)  % Compute iA
    [x_current,y_current,InfinityFlag_current] = PointAddition(x_current,y_current,InfinityFlag_current,x_A,y_A,InfinityFlag_A,a,p);
end
x_iA = x_current; y_iA = y_current; 

x_current = x_B; y_current = y_B; InfinityFlag_current = false;
for m = 1:(j-1)  % Compute jB
    [x_current,y_current,InfinityFlag_current] = PointAddition(x_current,y_current,InfinityFlag_current,x_B,y_B,InfinityFlag_B,a,p);
end
x_jB = x_current; y_jB = y_current; 

[x_final,y_final,~] = PointAddition(x_iA,y_iA,false,x_jB,y_jB,false,a,p)

% Problem 7
p = 27001; alpha = 101; a_u = 21768; a_v = 9898;
b_u = SquareAndMultiply(alpha,a_u,p)
b_v = SquareAndMultiply(alpha,a_v,p)
K = SquareAndMultiply(b_v,a_u,p)
K = SquareAndMultiply(b_u,a_v,p)
K = SquareAndMultiply(alpha,a_v*a_u,p)

function [x3,y3,InfinityFlag3] = PointAddition(x1,y1,InfinityFlag1,x2,y2,InfinityFlag2,a,p)
% Elliptic curve point addition
    if InfinityFlag1 == true && InfinityFlag2 == false
        x3 = x2;
        y3 = y2;
        InfinityFlag3 = false;
    elseif InfinityFlag1 == true && InfinityFlag2 == true
        InfinityFlag3 = true;
    elseif InfinityFlag1 == false && InfinityFlag2 == true
        x3 = x1;
        y3 = y1;
        InfinityFlag3 = false;
    elseif InfinityFlag1 == false && InfinityFlag2 == false
        if x2 == x1 && y2 == -y1
            InfinityFlag3 = true;
        else
            InfinityFlag3 = false;
            if x1 == x2
                if mod(y1+y2,p) == 0
                    InfinityFlag3 = true;
                    x3 = 0;
                    y3 = 0;
                elseif y1 == y2
                    lambda = mod((3*x1^2+a)*MultiplicativeInverse(2*y1,p),p);
                    x3 = mod(lambda^2-x1-x2,p);
                    y3 = mod(lambda*(x1-x3)-y1,p);
                else
                    lambda = mod((y2-y1)*MultiplicativeInverse(x2-x1,p),p);
                    x3 = mod(lambda^2-x1-x2,p);
                    y3 = mod(lambda*(x1-x3)-y1,p);
                end
            else
                lambda = mod((y2-y1)*MultiplicativeInverse(x2-x1,p),p);
                x3 = mod(lambda^2-x1-x2,p);
                y3 = mod(lambda*(x1-x3)-y1,p);
            end   
        end
    end
end

function a_inv = MultiplicativeInverse(a,b)
% Computes a_inv mod b using Extended Euclidean Algorithm
% s*a + t*b = r = gcd(a,b)
% Vector inputs possible
% Negative inputs possible
    
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