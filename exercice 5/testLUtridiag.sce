s=100

rand("seed",s)
n=5;

A=rand(n,n);

GB = zeros(n,n);
for k=1:n
    GB(k,k) = A(k,k) ;
    if k<n
        GB(k+1,k) = A(k+1,k);
        GB(k,k+1) = A(k+1,k);
    end
end

disp(GB)

[L,U]= FactoLU(GB);

/*
U=triu(B)
disp(U);
L=tril(B)
disp(L);
*/
