clear; clc;

n = 2:30;
m = floor(n/2);

nm = zeros(1,length(n));
for i = 1:length(n)
    nm(i) = nchoosek(n(i),m(i));
end

plot(n,nm);

