function [L,U]=FactoLU(A)

[m,n]=size(A);
if m~=n then, error('la matrice n''est pas carree'), end;
zero=1.e-12;
for k=1:n-1
   if abs(A(k,k))<zero then, error('erreur : pivot nul'), end;
   for i=k+1:n
     A(i,k)=A(i,k)/A(k,k);
     for j=k+1:n
        A(i,j)=A(i,j)-A(i,k)*A(k,j);
     end;
   end;
end;
U=triu(A);
disp(U);
L=A-U+diag(ones(n,1));
disp(L);
endfunction;



