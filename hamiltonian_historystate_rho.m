X = [0 1; 1 0];
Z = [1 0; 0 -1];
I = eye(2);

phi = pi/8;
D = sin(phi) * X + cos(phi) * Z;

U1 = X;
U2 = D;

x = 0;

Hin = (1/4) * kron(I - (-1)^x * Z, kron(I + Z, I));

Hclock = (1/4) * kron(I, kron(I + Z, I - Z));

Hclockinit = (1/2) * kron(I, kron(I - Z, I));

Hprop1 = (1/2) * kron(I, kron(I, I + Z)) - (1/2) * kron(U1, kron(X, I) + kron(X, Z)); 
Hprop2 = (1/2) * kron(I, kron(I - Z, I)) - (1/2) * kron(U2, kron(I, X) - kron(Z, X));

Hprop = Hprop1 + Hprop2;

Hout = (1/4) * kron(I + Z, kron(I, I - Z));

H = Hin + Hclock + Hprop;

[v,e] = eig(H);
disp(e(1));
state = v(1:8,1);

III = eye(8);
IIZ = kron(I, kron(I, Z));
IZI = kron(I, kron(Z, I));
IZZ = kron(I, kron(Z, Z));
XIX = kron(X, kron(I, X));
XXI = kron(X, kron(X, I));
XXZ = kron(X, kron(X, Z));
XZX = kron(X, kron(Z, X));
ZII = kron(Z, kron(I, I));
ZIX = kron(Z, kron(I, X));
ZIZ = kron(Z, kron(I, Z));
ZZI = kron(Z, kron(Z, I));
ZZX = kron(Z, kron(Z, X));

ops = cell(1, 11);
ops{1} = III;
ops{2} = ZII;
ops{3} = ZZI;
ops{4} = IZZ;
ops{5} = XXI;
ops{6} = XXZ;
ops{7} = XIX;
ops{8} = XZX;
ops{9} = ZIX;
ops{10} = ZZX;
ops{11} = ZIZ;
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
%dp = 0.05 : 0.0085 : 0.15;
results1 = zeros(1, numSamples);

tic;
for run = 1 : numSamples
    disp(run);
    tot = 0;

    rho = state * state';
    for i = 1 : 3
        rho = depolarizeX2(dp(run), 3, i, rho);
    end
    
    pacc = 0;
    for i = 1 : 11
        prob = p(i) * trace(projectors{i} * rho);
        pacc = pacc + prob;
    end

    results1(run) = pacc;
end
toc;
disp(results1);

% xq = 0 : 0.001 : dp(numSamples);
% vq2 = interp1(dp,results1,xq,'spline');
% plot(dp,results1,'o',xq,vq2,':.');
% xlim([0 dp(numSamples)]);
% title('Spline Interpolation');