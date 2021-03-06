%Khu dieu kien bien bai toan thanh thang trong khong gian 3 chieu
function [kk,ff]=khudkbien(kk,ff,bcdof,bcval)
n=length(bcdof);
sdof=size(kk);
for i=1:n
    c=bcdof(i);
    for j=1:sdof
        kk(c,j)=0;
        kk(j,c)=0;
    end 
    kk(c,c)=1;
    ff(c)=bcval(i);
end