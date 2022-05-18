clear; clc;

p = 97; a = 1; b = 12; 
x = 0:p-1;
status = zeros(size(x));
status_zero = zeros(size(x));
y_square = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    y_square(i) = mod(xi^3+a*xi+b,p);
    if SquareAndMultiply(y_square(i),(p-1)/2,p) == 1
        status(i) = 1;
    elseif SquareAndMultiply(y_square(i),(p-1)/2,p) == 0
        status_zero(i) = 1;
    end
end
Nb_pts = nnz(status)*2+nnz(status_zero)+1

x1 = 0; y1 = 20; InfinityFlag1 = false;
x_orig = x1; y_orig = y1; InfinityFlag_orig = InfinityFlag1;
answ = 53; % private key/answer of DLP
for i = 2:answ
    [x1,y1,InfinityFlag1] = PointAddition(x1,y1,InfinityFlag1,x_orig,y_orig,InfinityFlag_orig,a,p);
end
disp(['x1 = ', num2str(x1), ', y1 = ', num2str(y1)])

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