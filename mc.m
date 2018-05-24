function y = mc(x0,p0,N)
for i = 1:N
    y(i) = x0 + sqrt(p0) * normrnd(0,1);
end
end

