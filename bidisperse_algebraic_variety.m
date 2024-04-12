
% activation model

clear
close all
clc

%% Kohler parameters


A = (1.4*1e-3); % in micro meters

kappa_parameter = 1.28; %non dim

rd1 = 0.065;
rd2 = 0.049;

B1 = (kappa_parameter*rd1^3);
B2 = (kappa_parameter*rd2^3);



f1 = @(x) A./x - B1./(x.^(3));
f2 = @(x) A./x - B2./(x.^(3));

%%
XX = [0.001:1e-6:10];
figure
plot(XX,100*f1(XX))
hold on
plot(XX,100*f2(XX))
axis([0 10 0 0.2])

%%

A1 = 6*1e-10;
alpha2 = 0.2199;
w = 0.3*1e6;

%% 
r_crit1 = sqrt(3*B1/A);
r_crit2 = sqrt(3*B2/A);

s_crit1 = 0.666666666666*A/r_crit1;
s_crit2 = 0.666666666666*A/r_crit2;

r_eq1 = sqrt(B1/(A-(A1*w/alpha2)));
r_eq2 = sqrt(B2/(A-(A1*w/alpha2)));


%%

x = [0.01:1e-3:4];
y = [0.01:1e-3:1];

%AV = zeros(length(x),length(y));
separatrix = -1*ones(2,length(x));
counter = 1;
for i=1:length(x)
    for j=1:length(y)
        AV = f1(x(i)) - f2(y(j));
        if abs(AV)<1e-5
            separatrix(:,counter) = [x(i);y(j)];
            counter = counter + 1;
        end
    end
100*i/length(x)
end

Sforeq = zeros(1,size(separatrix,2));
for i=1:size(separatrix,2)
    Sforeq(i) = A1*w - f1(separatrix(1,i))*(alpha2*0.5*(separatrix(1,i)+separatrix(2,i)));
end

figure
plot(abs(squeeze(Sforeq)),'.')

[~,Index] = min(abs(squeeze(Sforeq)));


figure;
plot(separatrix(1,:),separatrix(2,:),'.');
hold on;
plot(y,y,'color','black');
plot(r_crit1,r_crit2,'.','Color','red','markersize',10);
plot(r_crit1*ones(1,length(y)),y,'Color','red');
plot(separatrix(1,Index),separatrix(2,Index),'.','Color','red','markersize',30);


coupled_req1 = separatrix(1,Index);
coupled_req2 = separatrix(2,Index);

S0_eq1 = f1(separatrix(1,Index));
S0_eq2 = f2(separatrix(2,Index));

