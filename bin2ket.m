function out = bin2ket(binvect, k0, k1)

out = 1;
for i = 1 : length(binvect)
    if (binvect(i) == 0)
        out = kron(out, k0);
    else
        out = kron(out, k1);
    end
end