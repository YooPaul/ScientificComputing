load('salmon.mat')

coeffs = polyfit(t, salmon, 1);
data = coeffs(1);
save('A1.dat','data','-ascii')
data = polyval(coeffs, 80);
save('A2.dat','data','-ascii')
salmon_predicted = polyval(coeffs,t);

rmsr = sqrt(sum((salmon_predicted - salmon).^2) / length(salmon));
save('A3.dat','rmsr','-ascii')

coeffs = polyfit(t, salmon, 2);
data = polyval(coeffs, 80);
save('A4.dat','data','-ascii')

coeffs = polyfit(t, salmon, 5);
data = polyval(coeffs, 80);
save('A5.dat','data','-ascii')

coeffs = polyfit(t, salmon, 20);
data = polyval(coeffs, 80);
save('A6.dat','data','-ascii')

% y = N0 * e^(r*t)
% ln(y) = r*t + ln(N0)

nlog_salmon = log(salmon);
coeffs = polyfit(t, nlog_salmon,1);
data = [exp(coeffs(2)), coeffs(1)];
save('A7.dat','data','-ascii')

data = polyval(coeffs, 80);
data = exp(data);
save('A8.dat','data','-ascii')

rms_error = @(v)(exp_model(v(1),v(2),v(3))); 
vmin = fminsearch(rms_error, [0.001,-0.01,10]);
a = vmin(1);
b = vmin(2);
c = vmin(3);

t1= 80;
data = exp(a*t1^2+b*t1+c);
save('A9.dat','data','-ascii')

yfine = spline(t, salmon, 80);
save('A10.dat','yfine','-ascii')