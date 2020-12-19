n = 50;
T0 = 5;
T1 = 20;


tol = 1e-8; maxit = 10000;

tic;
[xJ, relresJ, resvecJ, itJ] = jacobi(A,b,tol,maxit);
toc
tic;
[xGS, relresGS, resvecGS, itGS] = gaussSeidel(A,b,tol,maxit);
toc
figure;
xtitle("Evolution de la convergence de la méthode jacobi et gaussSeidel", "Nombre d itérations", "Résidu");

plot(log10(resvecJ), 'k.');
plot(log10(resvecGS), 'r.');
legend(["jacobi";"gauss seidel"]);

