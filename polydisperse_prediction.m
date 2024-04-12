
% activation model

clear
close all
clc

%% Kohler parameters


A = (1.4*1e-3); % in micro meters
%B = (3.5*1e-4);

kappa_parameter = 1.28; %non dim

RHeq = @(x,rd) 1 + A./x - (kappa_parameter*rd^3)./(x.^(3));


%% Critical radius and supersaturation

r2_crit = 3*(kappa_parameter*0.065^3)/A;

s_crit = RHeq(sqrt(r2_crit),0.065);

rm = 0.9*r2_crit;

%% model parameters (micrometers)


N = 1e-9; % 1e-9 per micon^3 = 1e3 per cm3

A1 = 6*1e-10;
A2 = 3.5*1e20;
A3 = 50;
rhow = 1e-15;

alpha2 = 4*pi*rhow*A2*A3*N;
D = A3;


%% Noise intensity and vertical velocity

%w = 0.34207*1e6; % micron per second
w = 0.34207*1e6; % micron per second
s = 0.0*w;

%% Mean B parameter

rd_mean = 0.065;
rd_var = 0.05;
B = kappa_parameter*rd_mean^3;

%% analytical distributions

gamma1 = 2*D*A1*w/(alpha2);
gamma2 = 2*D*A;
gamma3 = 2*D*B;
gamma4 = 2*D*A1*s/alpha2;

rho = @(x) x.*exp(2*( (2/3)*(gamma1-gamma2)*sqrt(x.^3)/(gamma4^2) + 2*(gamma3/(gamma4*gamma4))*x.^(1/2)     ));
%rhos = @(x) x.*exp(2*( (2/3)*(gamma1-gamma2)*sqrt(x.^3)/(gamma4^2) + 2*(gamma3/(gamma4*gamma4))*x.^(1/2)  - 0.25*log(x)   ));


%% Integration paramters

dt = 1e-2;
tmax = 5*1e2;
n = floor(tmax/dt);

%% Initialise mean field solution

n_mf = 1e7;

R = zeros(1,n_mf);
S = zeros(1,n_mf);
R(1) = B/A;

NOISE = randn(1,n_mf);
    
for i=1:length(R)-1
    S(i+1) = S(i) - dt*alpha2*S(i)*sqrt(R(i)) + dt*A1*w + sqrt(dt)*A1*s*NOISE(i);
    R(i+1) = max(0,R(i) + dt*2*A3*(S(i)+1-RHeq(sqrt(R(i)),rd_mean)));
end

%% Initialise population of particles

numparts = 5*1e2;
rd_pop = rd_mean + rd_var*rd_mean*randn(1,numparts);

SS = zeros(1,n);
RR = zeros(numparts,n);

RR(:,1) = (kappa_parameter*rd_pop.^3)/A;

for i=1:n-1
    SS(i+1) = SS(i) - dt*alpha2*SS(i)*mean(sqrt(RR(:,i)),1) + dt*A1*w + sqrt(dt)*A1*s*NOISE(i);
for K = 1:numparts
    RR(K,i+1) = max(0,RR(K,i) + dt*2*A3*(SS(i)+1-RHeq(sqrt(RR(K,i)),rd_pop(K))));
end
100*i/n
end

%%

mR = zeros(1,n);
for i=1:n

mR(i) = mean(RR(:,i));

end


figure
subplot(1,2,1)

[H,Bins] = hist(rd_pop,30);
plot(Bins,H,'linewidth',2,'color','black')
hold on
%plot((3*(kappa_parameter*rd_mean^3)/A)*ones(1,100),linspace(0,50,100))
set(gca, 'fontsize', 22)
axis square 
box on
ylabel('counts','Interpreter','latex')
xlabel(['Dry radius $(\mu \mathrm{m})$'],'Interpreter','latex')

subplot(1,2,2)
hold on
for i=1:min(100, numparts)
plot(dt*[1:n],RR(i,:),'LineWidth',1.5,'color',[0.5 0.5 0.5])
end
for i=1:min(100, numparts)
plot(dt*[1:n],(3*(kappa_parameter*rd_pop(i)^3)/A)*ones(1,n),'color',[1 0.7 0.7],'LineWidth',2)    
end
plot(dt*[1:n],(3*(kappa_parameter*rd_mean^3)/A)*ones(1,n),'color','red','LineWidth',2)
plot(dt*[1:n],R(1:n),'LineWidth',5,'color','black')
plot(dt*[1:n],mR,'LineWidth',5,'color','green','LineStyle','--')
axis([0 tmax 0 2.5])
set(gca, 'fontsize', 22)
axis square 
box on
ylabel('Squared radius $(\mu \mathrm{m}^2)$','Interpreter','latex')
xlabel('time $(s)$','Interpreter','latex')



%% plot


x = [0.2:1e-5:3];
sx = [0.0008:1e-6:0.0014];


rd_for_estimation = rd_mean + rd_var*rd_mean*randn(1,1e7);

critical_supersaturation_sample = 2*A./(3*sqrt(3*kappa_parameter*(rd_pop.^3)/A));

r0_equilibria = (kappa_parameter*(rd_for_estimation.^3))/(A - (A1*w/alpha2));
S0_equilibria = (A1*w/alpha2)*sqrt((r0_equilibria).^(-1));

mean_r0_equilibrium = (kappa_parameter*(rd_mean.^3))/(A - (A1*w/alpha2));
mean_S0_equilibrium = (A1*w/alpha2)*sqrt((mean_r0_equilibrium).^(-1));

[pdf_crits,bins_crits] = hist(3*kappa_parameter*(rd_for_estimation.^3)/A,x);
dbins_crits = bins_crits(2)-bins_crits(1);
normalization_factor_crits = dbins_crits*sum(pdf_crits);

[pdf_scrits,bins_scrits] = hist(2*A./(3*sqrt(3*kappa_parameter*(rd_for_estimation.^3)/A)),sx);
dbins_scrits = bins_scrits(2) - bins_scrits(1);
normalization_factor_scrits = dbins_scrits*sum(pdf_scrits);

[pdf_eqs,bins_eqs] = hist(r0_equilibria,x);
dbins_eqs = bins_eqs(2)-bins_eqs(1);
normalization_factor_eqs = dbins_eqs*sum(pdf_eqs);

[pdf_seqs,bins_seqs] = hist(S0_equilibria,100);
dbins_seqs = bins_seqs(2)-bins_seqs(1);
normalization_factor_seqs = dbins_seqs*sum(pdf_seqs);

[pdf_mf,bins_mf] = hist(R,100);
dbins_mf = bins_mf(2)-bins_mf(1);
normalization_factor_mf = dbins_mf*sum(pdf_mf);


analytical_pdf = rho(x);
analytical_normalization_factor = (x(2)-x(1))*sum(analytical_pdf);


figure
subplot(1,3,1)

plot(bins_crits,pdf_crits/normalization_factor_crits,'linewidth',1,'color','red')
hold on
plot(bins_eqs,pdf_eqs/normalization_factor_eqs,'linewidth',1,'color','magenta')
plot(bins_mf,pdf_mf/normalization_factor_mf,'linewidth',1,'color','blue')
plot(x,analytical_pdf/analytical_normalization_factor,'linestyle','--','linewidth',2,'color','black')
axis square
%xlim([0.2 1])
box on
set(gca, 'fontsize', 20)
xlabel('Crtical squared radii $(\mu \mathrm{m}^2)$','interpreter','latex')
ylabel('Normalized frequency','Interpreter','latex')

subplot(1,3,2)
plot(100*bins_seqs,pdf_seqs/normalization_factor_seqs,'color','magenta')
hold on
plot(100*bins_scrits,pdf_scrits/normalization_factor_scrits,'color','red')
box on
axis square
set(gca, 'fontsize', 20)
xlabel('Crtical supersaturation $(\%)$','interpreter','latex')
ylabel('Normalized frequency','Interpreter','latex')

subplot(1,3,3)
hold on
for i=1:min(numparts,200)
plot(dt*[1:n],100*(2*A./(3*sqrt(3*kappa_parameter*(rd_pop(i).^3)/A)))*ones(1,n),'color',[0.7,0.7,0.7])
end
plot(dt*[1:n],100*S(1:n),'color','blue','linewidth',2)
plot(dt*[1:n],100*SS,'color','black','LineWidth',2)

axis square
axis([dt dt*n 0.1 0.12])
box on
set(gca, 'fontsize', 20)
xlabel('time $(s)$','interpreter','latex')
ylabel('Supersaturation $\%$','Interpreter','latex')
%%

difference_pdfs = abs(pdf_crits/normalization_factor_crits - analytical_pdf/analytical_normalization_factor);

theoretical_supersaturation_prediction = length(find(100*2*A./(3*sqrt(3*kappa_parameter*(rd_pop.^3)/A))<100*mean_S0_equilibrium))/numparts
polydisperse_total_activation = length(find(RR(:,end)>0.75))/numparts
Indeces = zeros(1,length(analytical_pdf));


% integration_threshold = 0.1;
% for i=1:length(x)
%     if x(i)<integration_threshold & x(i+1)>integration_threshold
%         Index_intersection1 = i
%     end
% end
% 
% 
% [~,loc_mode1] = max(pdf_crits);
% [~,loc_mode2] = max(analytical_pdf);
% 
% 
% Frac1 = 0;
% Frac2 = 0;
% 
% for i=Index_intersection1:length(analytical_pdf)
% 
%     Frac1 = Frac1 + (x(2)-x(1))*(analytical_pdf(i)/analytical_normalization_factor);
% 
% end
% 
% 
% 
% theoretical_total_activation = Frac1
% polydisperse_total_activation = length(find(RR(:,end)>0.75))/numparts
% %apriori_total_activation = length(find(kappa_parameter*(rd_pop.^3)/A>kappa_parameter*(rd_mean^3)/A))/numparts
% 
% 
% 
% 
% 
% 
% 
% 
