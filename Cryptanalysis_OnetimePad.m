clear;

Ciphertext1 = 'bccxfcuebxmymiaaydkjazaytgujkbgwrbsgrgylzxknuouxndzlbxpfesrnaerxnoagbeseppdnntymrvwvnbwgjylnlvhipcsyyugqugzmyksrblwjrowiswbdonwodvcsntkgqyrtbho';
Ciphertext2 = 'jcecraxrbkfvaiieyfmanhtunnidbbeirctnnyamwxmduoujtxswrjogckrhalksxffukkrtrjlpmlqlgbmoiiupjeocifotugwsxleg';

L = min(length(Ciphertext1),length(Ciphertext2));
Diff = mod(Ciphertext1(1:L) - Ciphertext2(1:L),26);

P1 = 'fourscoreandsevenyearsagoourfathersbroughtforthonthiscontinentanewnationconceivedinLibertyanddedicatedtothepropositionthatallmenarecreatedequal';
P2 = 'nowweareengagedinagreatcivilwartestingwhetherthatnationoranynationsoconceivedandsodydicatedcanlongendure';
if length(P1)>=length(P2)
    L2 = min(length(P1),length(Diff));
    P2 = char(mod(P1(1:L2)-97 - Diff(1:L2),26)+97)
else
    L2 = min(length(P2),length(Diff));
    P1 = char(mod(P2(1:L2)-97 + Diff(1:L2),26)+97)
end

