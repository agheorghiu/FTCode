numPhysical = 5;

I = speye(2);
X = sparse([0 1; 1 0]);
Z = sparse([1 0; 0 -1]);

k0 = sparse([1; 0]);
k1 = sparse([0; 1]);

kp = sparse([1/sqrt(2); 1/sqrt(2)]);
km = sparse([1/sqrt(2); -1/sqrt(2)]);

mZ = sparse(zeros(2^numPhysical, 2^numPhysical));
mI = speye(2^numPhysical);

for i = 0 : 2^numPhysical - 1
    binvect = de2bi(i, numPhysical);
    kz = bin2ket(binvect, k0, k1);
    kx = bin2ket(binvect, kp, km);
    if (sum(binvect) > 1)
        mZ = mZ - kz * kz';
    else
        mZ = mZ + kz * kz';
    end
end

mX = 1;
for i = 1 : numPhysical
    mX = kron(mX, X);
end

ops = cell(1, 11);
ops{1} = kron(mI, kron(mI, mI));
ops{2} = kron(mZ, kron(mI, mI));
ops{3} = kron(mZ, kron(mZ, mI));
ops{4} = kron(mI, kron(mZ, mZ));
ops{5} = kron(mX, kron(mX, mI));
ops{6} = kron(mX, kron(mX, mZ));
ops{7} = kron(mX, kron(mI, mX));
ops{8} = kron(mX, kron(mZ, mX));
ops{9} = kron(mZ, kron(mI, mX));
ops{10} = kron(mZ, kron(mZ, mX));
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


numSamples = 101;
inc = 0.01;
dp = (0 : (numSamples-1)) * inc;
results2 = zeros(1, numSamples);
encstate = sparse(csvread('xencoded.csv'));

tic;
for run = 1 : numSamples
    disp(run);
    rho = encstate' * encstate;
    
    for i = 1 : 3*numPhysical
        rho = depolarizeX2(dp(run), 3*numPhysical, i, rho);
    end
    
    pacc = 0;
    for i = 1 : length(p)
        prob = p(i) * trace(projectors{i} * rho);
        pacc = pacc + prob;
    end
    
    results2(run) = pacc;
end
toc;
disp(results2);