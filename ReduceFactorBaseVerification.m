clear; clc; close all;

Factor_base = primes(100); cardinal_base = length(Factor_base); B = max(Factor_base);
N = 101;

% Reduce factor base
Factor_base_par = Factor_base; % Variable exists only inside parallel for-loop to minimize overhead
keep_in_factor_base = zeros(1,cardinal_base);
keep_in_factor_base(1) = 1; % Keep 2 in the factor base
for j = 2:cardinal_base
    m = mod(N,Factor_base_par(j));
    m_init = mod(N,Factor_base_par(j));
    for k = 2:(Factor_base_par(j)-1)/2
        m = mod(m*m_init,Factor_base_par(j));
    end
    if m == 1
        keep_in_factor_base(j) = 1;
    end
%     if mod(mod(N,Factor_base_par(j))^((Factor_base_par(j)-1)/2),Factor_base_par(j)) == 1
%         keep_in_factor_base(j) = 1;
%     end
end
New_factor_base_par = nonzeros(Factor_base_par.*keep_in_factor_base)';
New_cardinal_base = length(New_factor_base_par);