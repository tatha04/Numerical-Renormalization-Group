function [Sz,Q] = dekode(kode)
Q = rem(kode,100)-50;
Sz = (kode-rem(kode,100))/100-50;
end