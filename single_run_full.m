
% activation model

clear
close all
clc

%% Kohler parameters


A = (1.4*1e-3); % in micro meters
B = (3.5*1e-4);

k = 1.28; %non dim

RHeq = @(x,coeff) 1 + A./x - B./(x.^(3));


%% Critical radius and supersaturation

r2_crit = 3*B/A;

s_crit = RHeq(sqrt(r2_crit));

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
w = 0.32*1e6; % micron per second
s = 0.4*w;

%% Seq

Seq = A1*w/(alpha2*sqrt(rm));


disp(strcat('Critical S:',num2str(100*(s_crit-1))))
disp(strcat('Equilibrium S:',num2str(100*Seq)))
disp(strcat('Critical r2:',num2str(r2_crit)))
disp(strcat('Steady r2:',num2str(rm)))
disp(strcat('Threshold 2alpha2*A/3:',num2str(100*2*alpha2*A/3)))
disp(strcat('SS timescale:',num2str(100*A1*w)))

epsilon = A1*s/(alpha2*sqrt(rm));




%% Initialise solution

dt = 1e-2;
tmax = 1e4;
n = floor(tmax/dt);

R = zeros(1,n);
RR = zeros(1,n);
RRR = zeros(1,n);
S = zeros(1,n);
R(1) = B/A;
RR(1) = B/A;
RRR(1) = B/A;
    
for i=1:n-1
    S(i+1) = S(i) - dt*alpha2*S(i)*sqrt(R(i)) + dt*A1*w + sqrt(dt)*A1*s*randn();
    R(i+1) = max(0,R(i) + dt*2*A3*(S(i)+1-RHeq(sqrt(R(i)),k)));
    RR(i+1) = max(0,RR(i) + dt*2*A3*(A1*w/(alpha2*sqrt(RR(i)))+1-RHeq(sqrt(RR(i)),k))+sqrt(dt)*2*A3*(A1*s/(alpha2*sqrt(RR(i))))*randn()  );
    RRR(i+1) = max(0,RRR(i) + dt*2*A3*(Seq+1-RHeq(sqrt(RRR(i)),k))+sqrt(dt)*2*A3*(A1*s/(alpha2*sqrt(rm)))*randn()  );
end

%% plot

figure
plot(R,'LineWidth',2)
hold on
plot(RR,'LineWidth',2)
plot(RRR,'LineWidth',2)
plot(r2_crit*ones(size(R)),'black')
ylim([0,1.2])
