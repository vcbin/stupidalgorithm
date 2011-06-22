clc;clear all;

num_dims=2;
tau_1=1/sqrt(2*sqrt(num_dims));
tau_2=1/sqrt(2*num_dims);
x=linspace(0,15);
mean=0;
sigma=(sqrt(num_dims)+1)/(2*num_dims);
h=figure;
set(h,'name','Lognormal distribution pdf','Numbertitle','off');
y=lognpdf(x,mean,sigma);
plot(x,y);
%ylabel('Lognormal distribution pdf');
legend(strcat('\sigma=',num2str(sigma)));
h=figure;
set(h,'name','Lognormal distribution cdf','Numbertitle','off');
y=logncdf(x,mean,sigma);
plot(x,y);
ylabel('Probability');
legend(strcat('\sigma=',num2str(sigma)));

h=figure;
set(h,'name','Normal distribution pdf','Numbertitle','off');
x=linspace(-1,1);
y=normpdf(x,0,0.15);
plot(x,y,'b');
hold on;
y=normpdf(x,0,0.3);
plot(x,y,'g');
legend('\sigma=0.15','\sigma=0.3');
% prob=quad(lognpdf(x,mean,sigma),10,inf)