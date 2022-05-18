clear;

p = 11;
Squaresmodp = zeros(1,ceil(p/2));

for i = 1:ceil(p/2)
    Squaresmodp(i) = mod(i^2,p);
    QR = unique(sort(Squaresmodp));
end