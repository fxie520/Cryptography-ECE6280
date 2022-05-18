clear; clc;

n = 317940011; b = 77537081;
qm = Euclidean(b,n); m = length(qm);
c = zeros(1,m+1);
c(1) = 1; c(2) = qm(1);
d(1) = 0; d(2) = 1;
phi = zeros(1,m+1);
r = zeros(2,m+1);
for j = 3:m+1
    c(j) = qm(j-1)*c(j-1) + c(j-2);
    d(j) = qm(j-1)*d(j-1) + d(j-2);
    phi(j) = (d(j)*b-1)/c(j);
    if abs(phi(j)-round(phi(j)))<0.001  % Numerical instabilities
        coeff = [1 -(n-phi(j)+1) n];
        r = roots(coeff);
        if angle(r(1))<0.001 && angle(r(2))<0.001  % If roots are positive real values
            if abs(abs(r(1))-round(abs(r(1))))<0.001 && abs(abs(r(1))-round(abs(r(1))))<0.001 && abs(r(1))<n && abs(r(2))<n
                p = r(1)
                q = r(2)
            end
        end
    end
end

function q = Euclidean(A,B)
    r(1) = A;
    r(2) = B;
    m = 2;
    while r(m) ~= 0
        q(m) = floor(r(m-1)/r(m));
        r(m+1) = r(m-1) - q(m)*r(m);
        m = m + 1;
    end
    q = q(2:end);
end
