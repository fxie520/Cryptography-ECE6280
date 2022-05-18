clear; clc;

a = 1; p = 71;  % Parameters of the elliptic curve

% Examples
x = 4; y = 5;  % Our point

Order = OrderOfPoint(x,y,a,p);

function Order = OrderOfPoint(x,y,a,p)
% Point must not be the infinity point
    InfinityFlag = false;
    Order = 1;
    x_initial = x; y_initial = y;
    while InfinityFlag == false
        [x,y,InfinityFlag] = PointAddition(x,y,x_initial,y_initial,a,p);
        Order = Order + 1;
    end
end

function [x3,y3,InfinityFlag] = PointAddition(x1,y1,x2,y2,a,p)
% Points 1 and 2 must not be the infinity point
    InfinityFlag = false;
    if x1 == x2
        if mod(y1+y2,p) == 0
            InfinityFlag = true;
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