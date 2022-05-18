clear; clc;

% Case 1
% alpha = 5; beta = 8563; p = 28703;
% Case 2
alpha = 10; beta = 12611; p = 31153;

n = p-1;
F = factor(n);
factorbase = unique(F);
N = numel(factorbase);  % Number of distinct prime factors
count = zeros(1,N);  % Exponent of each distinct prime factor
for i = 1:N
   count(i) = sum(factorbase(i)==F);
end
modulus = factorbase.^count;

remainders = zeros(1,N);
for i = 1:N
    remainders(i) = PohligHellman(p,alpha,beta,factorbase(i),count(i));
end

% Final answer
log_alpha_beta = ChineseRemainderTheorem(remainders,modulus)



function x = ChineseRemainderTheorem(remainders,m)
% Compute mod(x,M)
% x = remainder(i) mod m(i)

    M = prod(m);
    Mi = M./m;
    y = MultiplicativeInverse(Mi,m);
    x = mod(sum(remainders.*Mi.*y),M);
end

function sum = PohligHellman(p,alpha,beta,q,c)
    n = p-1;
    
    j = 0;
    beta_array = zeros(1,c);  % Array containing beta_0 to beta_(c-1)
    beta_array(1) = beta;
    a = zeros(1,c);  % Array containing a_0 to a_(c-1)
    while j <= c-1
        % Compute delta
        delta = SquareAndMultiply(beta_array(j+1),n/q^(j+1),p);

        % Find i
        alpha_nq = SquareAndMultiply(alpha,n/q,p); % Compute mod(alpha^(n/q),p)
        i = 0;
        value = 1;  % value of mod(alpha^(i*n/1))
        equal_flag = false;
        while i <= q-1 && not(equal_flag)
            if delta == value
                equal_flag = true;
                a(j+1) = i;
            end
            i = i+1;
            value = mod(value*alpha_nq,p);
        end

        % Compute beta_(j+1)
        alpha_inv = MultiplicativeInverse(alpha,p);
        alpha_inv_exponent = SquareAndMultiply(alpha_inv,a(j+1)*q^j,p);
        beta_array(j+2) = mod(beta_array(j+1)*alpha_inv_exponent,p);

        j = j+1;
    end
    sum = 0;  % mod(a,q^c)
    for i = 1:c
        sum = sum + a(i)*q^(i-1);
    end
end

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