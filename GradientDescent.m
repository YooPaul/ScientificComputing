f = @(V)((1 - V(1)).^2 + 10*(V(2) - V(1).^2).^2);
data = f([-1,8]);
save('A1.dat','data','-ascii')

data = fminsearch(f,[-1;8]);
save('A2.dat','data','-ascii')

df_dx1 = @(V)(-2*(1-V(1))-40*V(1)*(V(2)-V(1).^2)); % good habit to use the .^ on even a single numeric element
df_dx2 = @(V)(20*(V(2)-V(1).^2));
gradient_f = [df_dx1([-1;8]);df_dx2([-1;8])];
save('A3.dat','gradient_f','-ascii')

data = norm(gradient_f, Inf);
save('A4.dat','data','-ascii')

steepest_descent  = @(t)([-1-t*df_dx1([-1;8]);8-t*df_dx2([-1;8])]);
data = steepest_descent(0.1);
save('A5.dat','data','-ascii')

data = f(data);
save('A6.dat','data','-ascii')

f_wrap = @(t)(f(steepest_descent(t)));
X = fminbnd(f_wrap,0,0.1);
save('A7.dat','X','-ascii')

data = steepest_descent(X);
save('A8.dat','data','-ascii')


% Make initial guess
% V =[t;x1;x2]
tole = 1e-4;
x0 = [-1;8];
count = 0;
gradient_f = @(V)([df_dx1([V(1);V(2)]);df_dx2([V(1);V(2)])]);
phi = @(t,V)([V(1);V(2)] - t * gradient_f([V(1);V(2)]));
f_wrap = @(t,V)(f(phi(t,V)));
while norm(gradient_f(x0), Inf) >= tole
       % Make the function phi(t) = x - t*grad f(x)
       % Find the minimum of f(phi(t)) between t=0 and t=0.1
       X = fminbnd(f_wrap, 0, 0.1,optimset,x0);
       % Make new guess = phi(t)
       x0 = phi(X,x0);
       count = count + 1;

end

save('A9.dat','x0','-ascii')
save('A10.dat','count','-ascii')
