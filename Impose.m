function E1 = Impose(E0,L,SMALL)
% Impose degeneracy over closely lying states 

E1 = E0; 
i = 1;
while (i<=L)
    j = i+1;
    while( (j <= L) && (abs(E0(j)-E0(i)) < SMALL))
        j=j+1;
    end
    E1(i:j-1) = mean(E0(i:j-1));
    i=j;
end

end