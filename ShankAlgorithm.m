clear; clc;

% Case 1
% p = 24691; alpha = 106; beta = 12375;

% Case 2
p = 458009; alpha = 6; beta = 248388;
n = p - 1;
m = ceil(sqrt(n));
list1 = zeros(m,2);
list2 = zeros(m,2);

% Compute mod(alpha^m,p)
alpham = 1;
for k = 1:m
    alpham = mod(alpham*alpha,p);
end

% Compute list 1
for j = 1:m
    list1(j,1) = j-1;
    value = 1; % list(j,2)
    for k = 1:j-1
        value = mod(value*alpham,p);
    end
    list1(j,2) = value;
end
list1_sorted = sortrows(list1,2);

alpha_inv = MultiplicativeInverse(alpha,p);

% Compute list 2
for i = 1:m
    list2(i,1) = i-1;
    value = mod(beta,p); % list2(i,2)
    for k = 1:i-1
        value = mod(value*alpha_inv,p);
    end
    list2(i,2) = value;
end
list2_sorted = sortrows(list2,2);

% Find collisions
[collision,index1,index2] = intersect(list1_sorted(:,2),list2_sorted(:,2),'stable');
j_final = list1_sorted(index1,1); i_final = list2_sorted(index2,1);
fprintf('Log_alpha^beta = %d \n',mod(m*j_final+i_final,n));

function a_inv = MultiplicativeInverse(a,b)
    % Extended Euclidean Algorithm: s*a + t*b = r = gcd(a,b)
    % a_inv: the inverse of a modulo b (if exists)
    
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
end
