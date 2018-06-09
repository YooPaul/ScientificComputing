function E = exp_model(a,b,c)
    load('salmon.mat')
    salmon_predicted = exp(a*t.^2+b*t+c);
    E = sqrt(sum((salmon_predicted - salmon).^2) / length(salmon));
end