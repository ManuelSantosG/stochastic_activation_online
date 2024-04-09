





close all


time_vals = 10.^(linspace(1,4,20));


figure
load("frac_activation_s1.mat")
p1 = semilogx(time_vals,100*frac_act,'LineWidth',2);
hold on
load("frac_activation_s2.mat")
p2 = semilogx(time_vals,100*frac_act,'LineWidth',2);
load("frac_activation_s3.mat")
p3 = semilogx(time_vals,100*frac_act,'LineWidth',2);
load("frac_activation_s4.mat")
p4 = semilogx(time_vals,100*frac_act,'LineWidth',2);

set(gca,'fontsize',20);
grid on
ylabel('Fraction of activation $(\%)$','Interpreter','latex')
xlabel('time $(s)$','interpreter','latex')
title(strcat('$w=0.196 (m/s)$, ','$(\mathrm{S}_{c}-\mathrm{S}_{eq})/\mathrm{S}_{c} = 0.8 \%$'),'Interpreter','latex')
legend([p1,p2,p3,p4],'$\sigma = 0.1w$ $(m/s)$','$\sigma = 0.2w$ $(m/s)$','$\sigma = 0.3w$ $(m/s)$','$\sigma = 0.4w$ $(m/s)$','location','northwest','interpreter','latex')
