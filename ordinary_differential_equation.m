
a = 5;
x_dot = @(x)(a * sin(x));
f = @(t)(2* atan(exp(a * t)/ (1 + sqrt(2))));


dt = 0.1;
T = 1;
t = 0:dt:T;
n = length(t);
y = zeros(1,n);
y(1) = pi / 4;

for k = 1:n-1
    y(k+1) = y(k) + dt * x_dot(y(k));
end

result = y(n);
save('A1.dat','result','-ascii')

max_error_1 = 0;
for k = 1:n
    error = norm(y(k) - f((k - 1) * dt), Inf);
    if error > max_error_1
        max_error_1 = error;
    end
end

save('A2.dat','max_error_1','-ascii')

dt = 0.01;
t = 0:dt:T;
n = length(t);
y = zeros(1,n);
y(1) = pi / 4;
for k = 1:n-1
    y(k+1) = y(k) + dt * x_dot(y(k));
end

result = y(n);
save('A3.dat','result','-ascii')

max_error_2 = 0;
for k = 1:n
    error = norm(y(k) - f((k - 1) * dt), Inf);
    if error > max_error_2
        max_error_2 = error;
    end
end

save('A4.dat','max_error_2','-ascii')

divide_error = max_error_1 / max_error_2;
save('A5.dat','divide_error','-ascii')

dt = 1;
T = 100;
t = 0:dt:T;
n = length(t);
y = zeros(1,n);
y(1) = pi / 4;

for k = 1:n-1
    y(k+1) = y(k) + dt * x_dot(y(k));
end

result = y(n);
save('A6.dat','result','-ascii')

% v [xk+1 ,k] 
%f = @(k)(y(k+1) - a * sin(y(k+1)) - y(k));  
for k = 2:n
    y(k) = fzero(@(x)(x - a * sin(x) - y(k-1)),3);
end

result = y(n);
save('A7.dat','result','-ascii')

g = -9.8;
L = 2;
mu = 0.5;

A = [0 1; 
     -mu g/L];
dt = 0.05;
T = 10;
t = 0:dt:T;
n = length(t);

x0 = [1;0];

f = @(t,y)([y(2); g/L*y(1) - mu * y(2)]);

y = zeros(2,n);
y(:,1) = x0;

for k = 1: n-1
   y(:,k+1) = y(:,k) + dt * f(t(k),y(:,k));     
end

result = y(1,n);
save('A8.dat','result','-ascii')

for k = 1: n-1
   y(:,k+1) = (eye(2) - dt*A) \ y(:,k);     
end

result = y(1,n);
save('A9.dat','result','-ascii')


[t_out, x_out] = ode45(f,t,x0);

result = x_out(n,1);
save('A10.dat','result','-ascii')