%% A) 
clear; clc; close all;

Factor_base = primes(100); cardinal_base = length(Factor_base); B = max(Factor_base);
Nb_digits_init = 8; Times_of_test = 10000; Nb_digits_final = 15; 
N = round(logspace(Nb_digits_init,Nb_digits_final,Times_of_test));

Runtime = zeros(1,Times_of_test);
Nb_smooth = zeros(1,Times_of_test); % Number of smooth numbers
Nb_smooth_th = zeros(1,Times_of_test); % Number of smooth numbers in theory (asymptotically)
Times_smooth_check = zeros(1,Times_of_test);
Division_count = 0;

parfor i = 1:Times_of_test
    Factor_base_par = Factor_base; % Variable exists only inside parallel for-loop to minimize overhead
    u = log(N(i))/log(B)/2;
    Times_smooth_check(i) = 1000;
%     Times_smooth_check(i) = round(cardinal_base*u^u);
    x = ceil(sqrt(N(i))):ceil(sqrt(N(i)))+Times_smooth_check(i)-1; % x values
    Nb_to_check = x.^2 - N(i); % Values to check
    Nb_to_check = nonzeros(Nb_to_check)';
    Is_smooth = zeros(1,length(Nb_to_check)); % Logical array containing results of smoothness
    tic;
    for j = 1:length(Nb_to_check)
        Index_factor_base = 1;
        while Index_factor_base <= cardinal_base
            while mod(Nb_to_check(j),Factor_base_par(Index_factor_base)) == 0
                Nb_to_check(j) = Nb_to_check(j)/Factor_base_par(Index_factor_base);
                Division_count = Division_count + 1;
            end
            Index_factor_base = Index_factor_base+1;
        end
        Is_smooth(j) = Nb_to_check(j)==1;
    end
    Nb_smooth(i) = nnz(Is_smooth);
    Nb_smooth_th(i) = Times_smooth_check(i)*u^(-u);
    Runtime(i) = toc;
end

semilogx(N,Runtime./Times_smooth_check);
xlim([min(N),max(N)]); ylim([0,max(Runtime./Times_smooth_check)]);
title('Runtime per smooth check vs N value','Fontsize',16);
xlabel('N','Fontsize',16); ylabel('Runtime per smooth check','Fontsize',16);
grid on;

Average_division_per_smooth_check = Division_count/sum(Times_smooth_check);
val = Average_division_per_smooth_check/cardinal_base

% Figure2 = figure(2);
% set(Figure2,'defaulttextinterpreter','latex');
% semilogx(N,Nb_smooth,N,Nb_smooth_th);
% xlim([min(N),max(N)]);
% title('\bf{Number of smooth values between $\sqrt{N}$ and $\sqrt{N}+14999$ vs N value}','Fontsize',16);
% set(Figure2,'defaulttextinterpreter','tex');
% xlabel('N','Fontsize',16); ylabel('Number of smooth values','Fontsize',16);
% grid on;

%% B)
clear; clc; close all;

Factor_base = primes(230); cardinal_base = length(Factor_base); B = max(Factor_base);
Nb_digits_init = 3; Times_of_test = 10000; Nb_digits_final = 13; 
N = round(logspace(Nb_digits_init,Nb_digits_final,Times_of_test)+rand*1000);
Is_even = not(mod(N,2));
N = N+Is_even; % Make all N odd (QR mod 2)

Runtime = zeros(1,Times_of_test);
Nb_smooth = zeros(1,Times_of_test); % Number of smooth numbers
Times_smooth_check = zeros(1,Times_of_test);
% Division_count = 0;

parfor i = 1:Times_of_test
    % Reduce factor base
    Factor_base_par = Factor_base; % Variable exists only inside parallel for-loop to minimize overhead
    keep_in_factor_base = zeros(1,cardinal_base);
    keep_in_factor_base(1) = 1; % Keep 2 in the factor base
    for j = 2:cardinal_base
        m = mod(N(i),Factor_base_par(j));
        m_init = mod(N(i),Factor_base_par(j));
        for k = 2:(Factor_base_par(j)-1)/2
            m = mod(m*m_init,Factor_base_par(j));
        end
        if m == 1
            keep_in_factor_base(j) = 1;
        end
    end
    New_factor_base_par = nonzeros(Factor_base_par.*keep_in_factor_base)';
    New_cardinal_base = length(New_factor_base_par);
    
    Times_smooth_check(i) = 5000;
    x = ceil(sqrt(N(i))):ceil(sqrt(N(i)))+Times_smooth_check(i)-1; % x values
    Nb_to_check = x.^2 - N(i); % Values to check
    Nb_to_check = nonzeros(Nb_to_check)';
    Is_smooth = zeros(1,length(Nb_to_check)); % Logical array containing results of smoothness
    Divisibility_table = zeros(length(Nb_to_check),New_cardinal_base);
    tic;
    for j = 1:length(Nb_to_check)
        for k = 1:New_cardinal_base
            if j<= New_factor_base_par(k)
                if mod(Nb_to_check(j),New_factor_base_par(k)) == 0
                    Divisibility_table(j,k) = 1;
                end
            elseif mod(j,New_factor_base_par(k)) == 0
                Divisibility_table(j,k) = Divisibility_table(New_factor_base_par(k),k);
            else
                Divisibility_table(j,k) = Divisibility_table(mod(j,New_factor_base_par(k)),k);
            end
            
            if Divisibility_table(j,k) == 1
                while mod(Nb_to_check(j),New_factor_base_par(k)) == 0
                    Nb_to_check(j) = Nb_to_check(j)/New_factor_base_par(k);
%                     Division_count = Division_count + 1;
                end
            end
        end
        Is_smooth(j) = Nb_to_check(j)==1;
    end
    Nb_smooth(i) = nnz(Is_smooth);
    Runtime(i) = toc;
end

% Average_division_per_smooth_check = Division_count/sum(Times_smooth_check);
% val = Average_division_per_smooth_check/log(log(B));

semilogx(log10(N),smooth(Runtime./Times_smooth_check));
xlim([min(log10(N)),max(log10(N))]); ylim([0,max(Runtime./Times_smooth_check)]);
title('Runtime per smooth check vs log(log(N))','Fontsize',16);
xlabel('Number of digits of N','Fontsize',16); ylabel('Runtime per smooth check','Fontsize',16);
grid on;