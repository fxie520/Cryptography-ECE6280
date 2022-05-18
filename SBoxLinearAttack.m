clear;

% Build look up table
table = zeros(16,8);
for i = 0:15
    table(i+1,1:4) = de2bi(i,4,'left-msb');
end
table(1,5:8) = de2bi(12,4,'left-msb');
table(2,5:8) = de2bi(13,4,'left-msb');
table(3,5:8) = de2bi(3,4,'left-msb');
table(4,5:8) = de2bi(7,4,'left-msb');
table(5,5:8) = de2bi(2,4,'left-msb');
table(6,5:8) = de2bi(4,4,'left-msb');
table(7,5:8) = de2bi(8,4,'left-msb');
table(8,5:8) = de2bi(6,4,'left-msb');
table(9,5:8) = de2bi(14,4,'left-msb');
table(10,5:8) = de2bi(1,4,'left-msb');
table(11,5:8) = de2bi(0,4,'left-msb');
table(12,5:8) = de2bi(15,4,'left-msb');
table(13,5:8) = de2bi(11,4,'left-msb');
table(14,5:8) = de2bi(5,4,'left-msb');
table(15,5:8) = de2bi(9,4,'left-msb');
table(16,5:8) = de2bi(10,4,'left-msb');

% Build linear approximation table
Freq_table = zeros(16,16);
for i = 1:16
    for j = 1:16
        if (i==1) && (j==1)
            Freq_table(i,j) = 16;
        else
            Columns = [de2bi(i-1,4,'left-msb') de2bi(j-1,4,'left-msb')];
            subtable = table(:,Columns~=0);
            Freq_table(i,j) = 16-nnz(mod(sum(subtable,2),2));
        end
    end
end