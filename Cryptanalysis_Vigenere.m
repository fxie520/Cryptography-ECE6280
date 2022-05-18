clear;

Ciphertext = 'GNWSCGCFLQAOQNMJRHDNGNGTLYTREFASIZSJYHLPKAJZYDLWSMMGKEQOEUPYFXLZSKLBOGJEFECVTRLLAMSPLVIQRAAVLWUSRDWJVRHLBWDLEYDLVZSYHWEMSVHCRATZSGDEWSCANDFEMRDWYVDUPQKWDDZWIQRBWLLGCJEGQPGMWSGKPOLAWDHGMLEGLJLOJPRZEUEMKAYTLYTKPSFIRYUATKZSLWDCGFSXCEWNWLEWNWDUWRHTLLHHNGLYVPCCIQRRGDHDRJOBTROIWSMMTZLPKEHVGFGWZBASVZJNEWSCMNLZLSNGOGNIGPCXFHNRKBBYCYOWTYLIRYZGTKAYJTLPQVESCCUAWPBOAUMSLOQPMXTKPKOOXWBEANPUSRULRZEUEFSNOPRLHHYYLIRYQMRYTTWAQORZEREFWRZZSDDDNAWPWHYJRDEFWRWSYFLHEGLPHCGKHDYBLHHHYJCDXC';

% % % % % Guess m value using the Kasiski test

% Segment length (>=3)
SegLength = 3;

% Number of segments
NbSeg = length(Ciphertext)-SegLength+1;

% All segments with uppercase letters automatically transformed into integers
Seg = zeros(NbSeg,SegLength);
for i = 1:NbSeg
    Seg(i,:) = Ciphertext(i:i+SegLength-1);
end

% Unique segments (repeated segments are removed)
UniqueSeg = unique(Seg,'Stable','rows');

% Number of unique segments
UniqueSegCount = zeros(size(UniqueSeg,1),1);
for i = 1:size(UniqueSeg,1)
    Diff = Seg-ones(size(Seg,1),1)*UniqueSeg(i,:);
    UniqueSegCount(i) = nnz(all(Diff==0,2));
end

% Repititive segments (those appear more than once)
RepititiveSeg = UniqueSeg(UniqueSegCount-1>0,:);

% Times each one of them appears
RepititiveSegCount = UniqueSegCount(UniqueSegCount~=1);

% Positions they appear
Position = zeros(size(RepititiveSeg,1),max(UniqueSegCount));
for i = 1:size(RepititiveSeg,1)
    Diff = Seg-ones(size(Seg,1),1)*RepititiveSeg(i,:);
    Position(i,1:RepititiveSegCount(i)) = find((all(Diff==0,2))>0)';
end

% Difference in position (position2 - position1, position3 - position2 etc.)
DiffPos = zeros(size(Position) - [0 1]);
for j = 1:size(DiffPos,2)
    DiffPos(:,j) = Position(:,j+1) - Position(:,j);
end
DiffPos(DiffPos<0) = 0;

% Suspected m values
m_Kasiski = zeros(size(DiffPos,1),1);
if j == 1
    m_Kasiski = DiffPos;
elseif j>1
    previousgcd = DiffPos(:,1);
    for i = 1:j-1
        previousgcd = gcd(previousgcd,DiffPos(:,i+1));
    end
    m_Kasiski = previousgcd;
end

% % % % % Guess m value using index of coincidence

% Times each letter appears
ni = zeros(26,1);
for i = 1:26
    ni(i) = numel(find(Ciphertext == char(64+i)));
end
n = length(Ciphertext);
kobs = dot(ni,ni-1)/n/(n-1);
kuniform = 1/26;
keng = 0.065;

% Suspected m value
m_IC = (keng-kuniform)/(kobs-kuniform);

% Check which value of m is correct

% Begin by trying m = 5 (Kasiski)
m = 5;

Substrings = reshape(Ciphertext,m,length(Ciphertext)/m);

% Times each letter appears in each substring
ni_2 = zeros(26,5);
for i = 1:5
    Substring = Substrings(i,:);
    for j = 1:26
        ni_2(j,i) = numel(find(Substring == char(64+j)));
    end
    pi_obs = ni_2/size(Substrings,2); % Probability observed
end

% Probability in English
pi = [0.082;0.015;0.028;0.043;0.127;0.022;0.020;0.061;0.070;0.002;0.008;0.040;0.024;0.067;0.075;0.019;0.001;0.060;0.063;0.091;0.028;0.010;0.023;0.001;0.020;0.001];

% Correlation factor, which is maximum and approximately 0.065 when we have
% the right key
Correlation = zeros(26,5);
for Shift = 0:25
    pi_shifted = [pi_obs(1+Shift:end,:);pi_obs(1:Shift,:)];
    for i = 1:5
        Correlation(Shift+1,i) = dot(pi,pi_shifted(:,i));
    end
end

for i = 1:5
    plot(0:25,Correlation(:,i));
    hold on;
end
plot(0:25,ones(1,26)*keng,'LineStyle','--');
xlabel('Shift value','FontSize',16);
ylabel('Correlation factor','FontSize',16);
set(gca,'FontSize',14); % Axis fontsize
hold off;

% Keys
Key = zeros(5,1);
DecryptedSubstrings = zeros(size(Substrings));
for i = 1:5
    Key(i) = find(Correlation(:,i) == max(Correlation(:,i)))-1;
    DecryptedSubstrings(i,:) = mod(Substrings(i,:)-65-Key(i),26)+65;
end

DecryptedSubstrings = char(DecryptedSubstrings);
Plaintext = reshape(DecryptedSubstrings,1,length(Ciphertext))


