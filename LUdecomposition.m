clear all; close all; clc;

A = [65 -15 -20; -15 39 -18; -20 -18 63];
[L,U,P] = lu(A);
data = U * P * L;
save('A1.dat','data','-ascii')

b = [50; 0; 75];
data = zeros(1, 100);

for v = 1:100
    b(2) = v;
    y = L \ (P * b);
    x = U \ y;
    data(v) = x(2);
end

save('A2.dat','data','-ascii')

for v = 1:100
    b(2) =v;
    x = inv(A) * b;
    data(v) = x(1);
end

save('A3.dat','data','-ascii')


s = sqrt(2) / 2; 
A = [-s 1 0 0 0 0 0 0 0 s 0 0 0     
     -s 0 0 0 0 0 0 0 -1 -s 0 0 0     
     0 -1 1 0 0 0 0 0 0 0 0 0 0     
     0 0 0 0 0 0 0 0 0 0 -1 0 0     
     0 0 -1 s 0 0 0 0 0 0 0 -s 0     
     0 0 0 -s 0 0 0 0 0 0 0 -s -1  
     0 0 0 -s -1 0 0 0 0 0 0 0 0
     0 0 0 0 1 -1 0 0 0 0 0 0 0  
     0 0 0 0 0 0 0 0 0 0 0 0 1  
     0 0 0 0 0 1 -1 0 0 -s 0 s 0
     0 0 0 0 0 0 0 0 0 s 1 s 0
     0 0 0 0 0 0 1 -1 0 0 0 0 0
     0 0 0 0 0 0 0 0 1 0 0 0 0]; 
 
b = [0;0;0;0;0;0;0;0;4;0;10;0;5];

[L,U,P] = lu(A);
y = L \ (P * b);
x = U \ y;
save('A4.dat','y','-ascii')
save('A5.dat','x','-ascii')

data = A \ b;
save('A6.dat','data','-ascii')


while max(abs(x)) <= 30
    b(9) = b(9) + 0.01;
    y = L \ (P * b);
    x = U \ y;    
end

data = b(9); 
save('A7.dat','data','-ascii')


A = [1e-20 1;1 1];
data = cond(A);
save('A8.dat','data','-ascii')

L = [1 0;1e20 1];
U = [1e-20 1;0 1-1e20];
data = L * U;
save('A9.dat','data','-ascii')

L = [1 0;1e-20 1];
U = [1 1;0 1-1e-20];
data = L * U;
save('A10.dat','data','-ascii')













