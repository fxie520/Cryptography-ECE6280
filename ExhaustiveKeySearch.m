clear;

f_ans = [0 0 0 0 1 1 0 0 1 1 0 0 1 0 0 1 0 1 1 1 0 1 1 0 0 0 1 0 0 1 0];
for i = 0:4095
   i_binary = de2bi(i,12,'left-msb');
   x = zeros(1,31); x(1:3) = i_binary(1:3);
   y = zeros(1,31); y(1:4) = i_binary(4:7);
   z = zeros(1,31); z(1:5) = i_binary(8:12);
   f = zeros(1,31);
   for j = 1:31
       if j>3
           x(j) = mod(x(j-3)+x(j-2),2);
       end
       if j>4
           y(j) = mod(y(j-4)+y(j-1),2);
       end
       if j>5
           z(j) = mod(z(j-5)+z(j-3),2);
       end
       f(j) = mod(x(j)*y(j)+y(j)*z(j)+z(j),2);
   end
   if nnz(f == f_ans) == 31
       x_ans = x;
       y_ans = y;
       z_ans = z;
   end
end
