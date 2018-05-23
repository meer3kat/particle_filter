function [x2f,w2f] = smc(x1,p1,w1,N,zi,m1,ml,mt,sigma1,sigmal,sigmat)

%mc_vector = normrnd(0,1,[1,N]);
%x2 = x1 + sqrt(p1) .* mc_vector;
x2 = x1;
x2f = x1;
n=@(x,m,s) (1./(sqrt(2*pi).*s)) .* exp(-((x-m).*(x-m))./(2.*s.*s));


w2 = w1 .* (n(zi,ml,sigmal) .* n(x1,mt,sigmat) ./n(x1,m1,sigma1));
w2 = w2./sum(w2);
w2f = w2;
%resampling

tot = 0;

for j = 1:1:N
    c(j) = tot + w2(j);
    tot = c(j);
end
i=1;
for j = 1:1:N
    u(j) = (1/N)*(rand+j-1);

    while u(j)>c(i)
        i = i+1;
    end
    x2f(j) = x2(i);
    w2f(j) = 1/N;
end



end
