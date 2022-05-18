clear;

Factor_base = primes(13); cardinal_base = length(Factor_base); B = max(Factor_base);
N = 100; % Numbers to factor

x = 11:50; 
Is_smooth = zeros(1,length(x)); % Logical array containing results of smoothness
for i = 1:length(x)
    Index_factor_base = 1;
    while Index_factor_base <= cardinal_base
        while mod(x(i),Factor_base(Index_factor_base)) == 0
            x(i) = x(i)/Factor_base(Index_factor_base);
        end
        Index_factor_base = Index_factor_base+1;
    end
    Is_smooth(i) = x(i)==1; 
end
Nb_smooth = nnz(Is_smooth);