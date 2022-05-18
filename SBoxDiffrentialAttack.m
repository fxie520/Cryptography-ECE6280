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

Diff_freq_table = zeros(16);
for i = 0:15 % i: Delta X
    Truth_table = zeros(16);
    Truth_table(:,1:8) = table;
    Truth_table(:,9:12) = mod(de2bi(i,4,'left-msb')-Truth_table(:,1:4),2);
    for j = 1:16
        for k = 1:16
            if Truth_table(j,9:12) == table(k,1:4)
                Truth_table(j,13:16) = table(k,5:8);
            end
        end
    end
    DeltaY = bi2de(mod(Truth_table(:,5:8)+Truth_table(:,13:16),2),'left-msb');
    for j = 0:15
        Diff_freq_table(i+1,j+1) = nnz(find(DeltaY==j));
    end
end