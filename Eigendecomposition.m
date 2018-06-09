clear all ;close all; clc;

P(1:100) = 2; 
main_diag = diag(P);
above_below_val(1:99) = -1;
above_diag = diag(above_below_val, 1);
below_diag = diag(above_below_val, -1);
A = main_diag + above_diag + below_diag;

p = zeros(100,1);
for j = 1:100
    p(j) = 2 * (1 - cos(56 * pi / 101)) * sin(56 * pi * j / 101);
end

T = A - main_diag;
M = -P\T;
lambdas = abs(eig(M));
data = max(lambdas);
save('A1.dat','data','-ascii')

x0 = ones(100, 1);
tole = 1e-4;
x1 = P\(-T*x0 + p);
count = 1;
while norm(x1-x0,Inf) >= tole
    x0 = x1;
    x1 = P\(-T*x0 + p);
    count = count + 1;
end

save('A2.dat','count','-ascii')

true_x = A \ p;
data = norm(x1 - true_x, Inf);
save('A3.dat','data','-ascii')

P = tril(A);
T = A - P;
M = -P\T;
lambdas = abs(eig(M));
gauss_eigen = max(lambdas);
save('A4.dat','gauss_eigen','-ascii')

x0 = ones(100, 1);
x1 = P\(-T*x0 + p);
count = 1;
while norm(x1-x0,Inf) >= tole
    x0 = x1;
    x1 = P\(-T*x0 + p);
    count = count + 1;
end

save('A5.dat','count','-ascii')
data = norm(x1 - true_x, Inf);
save('A6.dat','data','-ascii')



w = 1.5
lower = tril(A);
upper = triu(A);
P = 1/w * main_diag + lower;
T = (w-1)/w * main_diag + upper;
M = -P\T;
lambdas = abs(eig(M));
data = max(lambdas);
save('A7.dat','data','-ascii')

minW = 1;
min_eigen = gauss_eigen
for w = 1.01:0.01:1.99
    P = 1/w * main_diag + lower;
    T = (w-1)/w * main_diag + upper;
    M = -P\T;
    lambdas = abs(eig(M));
    data = max(lambdas);
    if(data < min_eigen)
        min_eigen = data;
        minW = w;
    end
end

save('A8.dat','minW','-ascii')
P = 1/w * main_diag + lower;
T = (w-1)/w * main_diag + upper;

x0 = ones(100, 1);
x1 = P\(-T*x0 + p);
count = 1;
while norm(x1-x0,Inf) >= tole
    x0 = x1;
    x1 = P\(-T*x0 + p);
    count = count + 1;
end
save('A9.dat','count','-ascii')
data = norm(x1 - true_x, Inf);
save('A10.dat','data','-ascii')