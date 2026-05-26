clear all
close all
clc


mu = 0;
sigma_vet = 1:3;
figure
for i=1:length(sigma_vet)
    sigma=sigma_vet(i);
x = linspace(-5*sigma,5*sigma,1000);
pdf = 1/(sqrt(2*pi)*sigma) * exp(-(x-mu).^2/(2*sigma^2));
plot(x, pdf, 'LineWidth', 2)
hold on
grid on
end
legend('sigma=1','sigma=2','sigma=3','Fontsize',14)
title('PDF Gaussiana (formula esplicita)')



% 
 aa=randn(100000,1);
 figure
 hist(aa,20)