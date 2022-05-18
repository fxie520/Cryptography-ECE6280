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

Diff_freq_table = zeros(8);
for i = 0:7 % i: Delta X
    Truth_table = zeros(8,12);
    Truth_table(:,1:6) = table;
    Truth_table(:,7:9) = mod(de2bi(i,3,'left-msb')-Truth_table(:,1:3),2);
    for j = 1:8
        for k = 1:8
            if Truth_table(j,7:9) == table(k,1:3)
                Truth_table(j,10:12) = table(k,4:6);
            end
        end
    end
    DeltaY = bi2de(mod(Truth_table(:,4:6)+Truth_table(:,10:12),2),'left-msb');
    for j = 0:7
        Diff_freq_table(i+1,j+1) = nnz(find(DeltaY==j));
    end
end