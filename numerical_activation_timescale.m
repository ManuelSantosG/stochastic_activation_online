
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

rm = B/A;

%% model parameters (micrometers)


N = 1e-9; % 1e-9 per micon^3 = 1e3 per cm3

A1 = 6*1e-10;
A2 = 3.5*1e20;
A3 = 50;
rhow = 1e-15;

alpha2 = 4*pi*rhow*A2*A3*N;
D = A3;


%% Noise intensity and vertical velocity

w = 0.196*1e6; % micron per second


%% Seq

Seq = A1*w/(alpha2*sqrt(rm));


disp(strcat('Critical S:',num2str(100*(s_crit-1))))
disp(strcat('Equilibrium S:',num2str(100*Seq)))
disp(strcat('Critical r2:',num2str(r2_crit)))
disp(strcat('Steady r2:',num2str(rm)))



%%

s_vals = [0.01:0.01:0.4];


epsilon_vals = A1*s_vals*w/(alpha2*sqrt(rm));

Num_ens = 100;

%% Initialise solution
dt = 1e-2;
tmax = 1e5;
n = floor(tmax/dt);

act_time = zeros(length(s_vals),Num_ens);

for K = 1:length(epsilon_vals)


epsilon = epsilon_vals(K);

for nens = 1:Num_ens
    R = zeros(1,n);
    R(1) = rm;
    
    for i=1:n-1
        
        R(i+1) = max(0,R(i) + dt*2*A3*(Seq+1-RHeq(sqrt(R(i)),k))+sqrt(dt)*A3*epsilon*randn()  );
        if R(i+1)>r2_crit
            break
        end

    end

    act_time(K,nens) = i*dt;
    nens/Num_ens
end

K/4
end

%%

save('numerical_activation_times.mat','act_time')
