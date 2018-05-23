%for particle filter

clear 
tic
%load data
load('SP500');
%take log price

z = log(data); 
N = 100; 
x = zeros(length(data),N);
p = zeros(length(data),N);
w = zeros(length(data),N);
%para = [u,k,th,lam,ro]
%       u, kappa, theta, lamda, ro
%model parameters
u = 0.2016;
k = 2.1924;
th = 0.0133;
lam = 0.294;
ro = -0.6143;
%delta t = number of trading days
dt = 1/252;
%for monte carlo simulation numbers



x0 = 0.0233;
p0 = 0.00001;
x0 = mc(x0,p0,N);
w0 = 1/N * ones(1,N);
m1 = x0 + k*(th-x0)*dt + lam*ro*(z(2)-z(1) - (u-0.5*x0)*dt);
sigma1 = lam*sqrt(x0*(1-ro*ro)*dt);
ml = z(1) + (u - 0.5*x0)*dt;
mt = (1/(1+0.5*lam*ro*dt))*(x0 + k*(th-x0)*dt + 0.5*lam*ro*x0*dt);
sigmal=sqrt(x0) * dt;
sigmat=(1/(1+0.5*lam*ro*dt))*lam.*sqrt(x0*dt);

%for our first data sets. 
Fk = 1 - k*(dt) + 0.5*ro*lam*dt;
Ak = k*th*dt-lam*ro*u*dt;
Qk = lam*lam*(1-ro*ro)*x0*dt;
Hk = -0.5*dt;
Bk = z(1) + u*dt;
Rk = x0 * dt;

x(1,:) = Fk*x0 + Ak; %first prediction of x
%zp(1) = Hk*xp(1) + Bk; %predicted z (log price)
p(1,:) = Fk*p0*Fk'+Qk;
w(1,:) = w0;
xf(1) = sum(x(1,:).*w(1,:));
logprice(1) = Hk *xf(1) + Bk;

[x(2,:),w(2,:)] = smc(x(1,:),p(1,:),w(1,:),N,z(1),m1,ml,mt,sigma1,sigmal,sigmat);


for i = 2:1:length(data)
    i;
    
    Fk = 1 - k*(dt) + 0.5*ro*lam*(dt);
    Ak = k*th*dt - lam*ro*u*dt + lam*ro*(z(i)-z(i-1));
    Qk = lam * lam * (1-ro*ro).*x(i-1,:)*dt;

    
    m1 = x(i-1,:) + k.*(th-x(i-1,:))*dt + lam*ro*(z(i)-z(i-1) - (u-0.5.*x(i-1,:))*dt);
    sigma1 = lam*sqrt(x(i-1,:).*(1-ro*ro).*dt);
    ml = z(i-1) + (u- 0.5*x(i-1,:))*dt;
    mt = (1/(1+0.5*lam*ro*dt)).*(x(i-1,:) + k*(th-x(i-1,:))*dt + 0.5*lam*ro*x(i-1,:)*dt);
    sigmal=sqrt(x(i-1,:)) * dt;
    sigmat=(1/(1+0.5*lam*ro*dt))*lam.*sqrt(x(i-1,:)*dt);
    x(i,:) = Fk * x(i-1,:) + Ak;
    p(i,:) = Fk * p(i-1,:) * Fk' + Qk;
    
    [x(i,:),w(i,:)] = smc(x(i,:),p(i,:),w(i-1,:),N,z(i),m1,ml,mt,sigma1,sigmal,sigmat);
    xf(i) = sum(x(i,:).*w(i,:));
    
    Hk = -0.5*dt;
    Bk = z(i-1) + u*dt;
    logprice(i) = Hk *xf(i) + Bk;
    

    
end
price = exp(logprice);
%{
f1 = figure('position', [0, 0, 700, 500]);
plot(dtime, z, 'b', dtime, zk, 'r-');
datetick('x','yyyy');
legend({'Prediction';'SP500'},'FontSize',20)
datetick('x','yyyy');
xlabel('Time','FontSize',20);
ylabel('Index','FontSize',20);
hold off
%}
save('pfdata100.mat')
f1 = figure('position', [0, 0, 700, 500]);
plot(dtime, xf);
set(gca, 'FontSize', 15)
datetick('x','yyyy');
legend({'volatility PF'},'FontSize',15)
datetick('x','yyyy');
xlabel('Time','FontSize',15);
ylabel('Volatility','FontSize',15);
hold off
%}

f2 = figure('position', [0, 0, 700, 500]);
plot(dtime, data,'b-*');
hold on
plot(dtime,price,'rd')
set(gca, 'FontSize', 15)
datetick('x','yyyy');
legend({'price';'model price'},'FontSize',15)
datetick('x','yyyy');
xlabel('Time','FontSize',15);
ylabel('price','FontSize',15);
hold off


%saveas(f1,'ekfvol.png');




%{
f2 = figure
histu = hist(res,[-0.1:0.01:0.1]);
x= [-0.1:0.01:0.1];
bar(x,histu);
xlabel('residual','FontSize',14);
ylabel('frequency','FontSize',14);
title('plot of residual for EKF','FontSize',18)
%saveas(f2,'res.png');

sume = sum(abs(res));
avere = mean(sume)/length(res);
finale = exp(avere);



finale =

    1.0091

%}
toc