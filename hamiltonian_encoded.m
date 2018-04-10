I = speye(2);
X = sparse([0 1; 1 0]);
Z = sparse([1 0; 0 -1]);
H = sparse([1 1; 1 -1]);

k0 = sparse([1; 0]);
k1 = sparse([0; 1]);

mZ = sparse(zeros(2^7, 2^7));
mI = speye(2^7);

for i = 0 : 2^7 - 1
    binvect = de2bi(i, 7);
    kz = bin2ket(binvect, k0, k1);
    
    % check syndromes
    syn1 = mod(sum([binvect(1) binvect(3) binvect(5) binvect(7)]), 2);
    syn2 = mod(sum([binvect(2) binvect(3) binvect(6) binvect(7)]), 2);
    syn3 = mod(sum([binvect(4) binvect(5) binvect(6) binvect(7)]), 2);
    
    % correction
    err = bi2de([syn1 syn2 syn3]);
    if (err > 0)
        binvect(err) = 1 - binvect(err);
    end
    
    % recovery
    res = mod(sum(binvect), 2);
    
    mZ = mZ + ((-1)^res) * (kz * kz');
end

ops = cell(1, 11);
ops{1} = kron(mI, kron(mI, mI));
ops{2} = kron(mZ, kron(mI, mI));
ops{3} = kron(mZ, kron(mZ, mI));
ops{4} = kron(mI, kron(mZ, mZ));
ops{5} = kron(mZ, kron(mZ, mI));
ops{6} = kron(mZ, kron(mZ, mZ));
ops{7} = kron(mZ, kron(mI, mZ));
ops{8} = kron(mZ, kron(mZ, mZ));
ops{9} = kron(mZ, kron(mI, mZ));
ops{10} = kron(mZ, kron(mZ, mZ));
ops{11} = kron(mZ, kron(mI, mZ));

III = ops{1};

x = 0;
phi = pi/8;
coeffs = [7/4 (1/4)*(1 - (-1)^x) -1/4*(-1)^x -1/4 -1/2 -1/2 (-1/2)*sin(phi) (1/2)*sin(phi) (-1/2)*cos(phi) (1/2)*cos(phi) -1/4];
p = abs(coeffs) / sum(abs(coeffs));

projectors = cell(1, 11);
for i = 1 : 11
    projectors{i} = (III - sign(coeffs(i)) * ops{i}) / 2;
end
clear ops;


numSamples = 12;
numRepeats = 2000;
%inc = 0.091;
%dp = (0 : (numSamples-1)) * inc;
dp = 0.05 : 0.0085 : 0.15;
results2 = zeros(1, numSamples);
stdevs2 = zeros(1, numSamples);
encstate = sparse(csvread('encodedstate.csv'));

tic;
spmd
    disp(labindex);
    pacc = zeros(1, numRepeats);
    for k = 1 : numRepeats
        ket = encstate;
        for i = 1 : 21
            ket = depolarizeK2(dp(labindex), 21, i, ket);
        end

        for i = 1 : length(p)
            nket = ket;
            if (i == 5 || i == 6)
                nket = applyHadamards((1 : 14), 21, nket);
            end
            if (i == 7 || i == 8)
                nket = applyHadamards((1 : 7), 21, nket);
                nket = applyHadamards((15 : 21), 21, nket);
            end
            if (i == 9 || i == 10)
                nket = applyHadamards((15 : 21), 21, nket);
            end
            prob = p(i) * (nket' * projectors{i} * nket);
            pacc(k) = pacc(k) + prob;
        end
    end
    tot = mean(pacc);
    sdev = std(pacc);
end
toc;
for i = 1 : 12
    results2(i) = tot{i};
    stdevs2(i) = sdev{i};
end
disp(results2);
disp(stdevs2);
save('data.mat', 'results2', 'stdevs2');