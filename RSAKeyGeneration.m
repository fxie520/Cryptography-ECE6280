clear; clf;

% Data points
Nb_bits = linspace(1024,20480,20)';
Average_Real_time = [0.168 0.766 2.709 5.241 9.859 18.642 26.405 80.119 67.788 136.290 114.363 319.469 302.065 307.313 928.456 1092.445 1192.511 1287.833 1670.106 2134.729]'/10;

% Linear least squares fit
phi1 = @(x) x.^1;
phi2 = @(x) x.^2;
phi3 = @(x) x.^3;

X = [phi1(Nb_bits), phi2(Nb_bits), phi3(Nb_bits)];
% X = [phi1(Nb_bits), phi2(Nb_bits)];
A = X' * X;
B = X' * Average_Real_time;

delta = A\B; % delta = (A^(-1)) * B also possible;

% Plot data points and best fit function
scatter(Nb_bits,Average_Real_time);
hold on;
x = linspace(0,20480,1000);
y = delta(1)*phi1(x) + delta(2)*phi2(x) + delta(3)*phi3(x);
% y = delta(1)*phi1(x) + delta(2)*phi2(x);
plot(x,y);
xlabel('Number of bits of p/q','fontsize',16); xlim([0 20480]);
ylabel('Average real (wall clock) time (s)','fontsize',16); 
grid on;

% Estimated wall clock time (32768 bits)
Estimation = delta(1)*32768+delta(2)*32768^2+delta(3)*32768^3
% Estimation = delta(1)*32768+delta(2)*32768^2
