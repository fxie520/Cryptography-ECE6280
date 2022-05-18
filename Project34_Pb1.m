clear; clc;

p = 10357; a = 39; b = 101;  % Parameters of the curve
x = 0:p-1; % Abscissa of points
y_square = zeros(size(x));
Is_QR = zeros(size(x));  % Whether y^2 is QR
Is_Zero = zeros(size(x));  % Whether y^2 is 0

for i = 1:length(x)  % Exhaustive search
    xi = x(i);
    y_square(i) = mod(xi^3+a*xi+b,p);
    if SquareAndMultiply(y_square(i),(p-1)/2,p) == 1  % Euler's criterion
        Is_QR(i) = 1;
    elseif SquareAndMultiply(y_square(i),(p-1)/2,p) == 0  % y^2 = 0 mod p
        Is_Zero(i) = 1;
    end
end
Nb_pts = nnz(Is_QR)*2+nnz(Is_Zero)+1

x_A = 117; y_A = 271;  % Point A

Order = OrderOfPoint(x_A,y_A,a,p)

x_B = 1651; y_B = 6391;  % Point B

% Method 1: Pollard-Rho

% Method 2 (for Verification): exhaustive point addition
m = 1;
x_initial = x_A; y_initial = y_A;
x_current = x_A; y_current = y_A;
while x_current~=x_B || y_current~=y_B
    [x_current,y_current,InfinityFlag] = PointAddition(x_current,y_current,x_initial,y_initial,a,p);
    m = m + 1;
end
disp("Method 2: m=");
disp(m);

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