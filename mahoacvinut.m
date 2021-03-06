function[index]=mahoacvinut(iel,node,nnel,ndof)
nd(1)=node(iel,1);
nd(2)=node(iel,2);
k=0;
for i=1:nnel
    start=(nd(i)-1)*ndof;
    for j=1:ndof
        k=k+1;
        index(k)=start+j;
    end
end