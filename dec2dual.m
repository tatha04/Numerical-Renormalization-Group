function z = dec2dual(M,N)
    z = zeros(N,1);
    i=1;
    while(M > 1)
        z(i) = mod(M,2);
        M = (M-z(i))/2;
        i = i+1;
    end
    if (M > 0)
        z(i) = M;
    end
end