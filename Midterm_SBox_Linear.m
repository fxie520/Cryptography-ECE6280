clear;

% Build look up table
table = zeros(8,6);
for i = 0:7
    table(i+1,1:3) = de2bi(i,3,'left-msb');
end
table(1,4:6) = de2bi(1,3,'left-msb');
table(2,4:6) = de2bi(7,3,'left-msb');
table(3,4:6) = de2bi(3,3,'left-msb');
table(4,4:6) = de2bi(4,3,'left-msb');
table(5,4:6) = de2bi(2,3,'left-msb');
table(6,4:6) = de2bi(0,3,'left-msb');
table(7,4:6) = de2bi(5,3,'left-msb');
table(8,4:6) = de2bi(6,3,'left-msb');

% Build linear approximation table
Freq_table = zeros(8);
for i = 1:8
    for j = 1:8
        if (i==1) && (j==1)
            Freq_table(i,j) = 8;
        else
            Columns = [de2bi(i-1,3,'left-msb') de2bi(j-1,3,'left-msb')];
            subtable = table(:,Columns~=0);
            Freq_table(i,j) = 8-nnz(mod(sum(subtable,2),2));
        end
    end
end