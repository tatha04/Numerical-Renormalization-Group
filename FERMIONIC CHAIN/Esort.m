function [Es,cs,Ks] = Esort(E0,c0,Kodep0,L)
[Es,ix] = sort(E0);
cs = zeros(L); Ks = zeros(L,1);
for i = 1:L
    cs(:,i) = c0(:,ix(i));
    Ks(i)   = Kodep0(ix(i));
end
end