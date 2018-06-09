
dt = 0.5;
T = 100;
t = 0:dt:T;
n = length(t);

a = 0.7;
b = 0.8;
tor = 12.5;

I = @(t)(1/10 * (5+sin(pi*t/10)));
% v = [v; w]
%v_dot = @(t,y)(y(1) - 1/3 * y(1)^3 - y(2) + I(t));
%w_dot = @(t,y)((a + y(1) - b * y(2))/ tor);
f = @(t,v)([v(1) - 1/3 * v(1).^3 - v(2) + I(t);(a+v(1)-b*v(2))/tor]);
y0 = [1; 0];
y = zeros(2, n);
y(:, 1) = y0;

for k = 1:n - 1
    y(:,k+1) = y(:,k) + dt*f(t(k)+dt/2,y(:,k)+dt/2*f(t(k),y(:,k)));
    %y(1, k+1) = y(1,k) + dt * v_dot(t(k) + dt/2, [y(1,k) + dt/2 * v_dot(t(k),y(:,k)); y(2,k) + dt/2 * w_dot(t(k),y(:,k))]);
    %y(2, k+1) = y(2,k) + dt * w_dot(t(k) + dt/2, [y(1,k) + dt/2 * v_dot(t(k),y(:,k));y(2,k) + dt/2 * w_dot(t(k),y(:,k))]);
    %[t_out,y] = ode23(wrapper,t,y0);
end

result = y(1,n);
save('A1.dat','result','-ascii')

[local_max,max_ind] = findpeaks(y(1,:));
[local_min,min_ind] = findpeaks(-y(1,:));

%[max_val,max_ind] = max(y(1,:));
%[min_val,min_ind] = min(y(1,:));

%hold on
%plot(t,y(1,:),'k')
%plot(t(max_ind),local_max,'r.')
%plot(t(min_ind),-local_min,'b.')

amp = local_max(1) + local_min(1);
period = t(max_ind(3)) - t(max_ind(1));

%amp = max_val(1) - min_val(1);
%period = t(max_ind(3)) - t(max_ind(1));

save('A2.dat','amp','-ascii')
save('A3.dat','period','-ascii')


% fourth order Runge-Kutta

for k = 1:n - 1
    f1 = f(t(k),y(:,k));
    f2 = f(t(k)+dt/2,y(:,k)+dt/2*f1);
    f3 = f(t(k)+dt/2,y(:,k)+dt/2*f2);
    f4 = f(t(k)+dt,y(:,k)+ dt*f3);
    y(:,k+1) = y(:,k) + dt/6 * (f1+2*f2+2*f3+f4);
    
end

result = y(1,n);
save('A4.dat','result','-ascii')

[local_max,max_ind] = findpeaks(y(1,:));
[local_min,min_ind] = findpeaks(-y(1,:));

amp = local_max(1) + local_min(1);
period = t(max_ind(3)) - t(max_ind(1));

%[max_val,max_ind] = max(y(1,:));
%[min_val,min_ind] = min(y(1,:));

%amp = max_val(1) - min_val(1);
%period = t(max_ind(3)) - t(max_ind(1));

save('A5.dat','amp','-ascii')
save('A6.dat','period','-ascii')

% boundary value problem

x0 = 1;
xT = 2;
T = 6;

dt = 0.01;
t = 0:dt:T;
N = floor(T/dt);

v2 = -2 * ones(N - 1, 1);
v1 = ones(N - 2, 1);
A = diag(v2) + diag(v1,1)+diag(v1,-1);
A = (1/dt^2) * A;
A = A + eye(N -1);

b = zeros(N - 1, 1);
for i = 1: N - 1
    b(i) = 4*cos(5*t(i+1));
end
b(1) = -x0 /dt^2 + 4 * cos(5*t(1));
b(end) = -xT/dt^2 + 4 * cos(5*t(end));

x_int = A\b;
x = [x0; x_int; xT];

int_points = length(x_int);
save('A7.dat','int_points','-ascii')

approx = x(floor(3/dt)+1);
save('A8.dat','approx','-ascii')

[max_val,max_time] = max(x);
[min_val,min_time] = min(x);

result = t(max_time(1));
save('A9.dat','result','-ascii')

result = t(min_time(1));
save('A10.dat','result','-ascii')
