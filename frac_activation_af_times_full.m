
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

w = 0.27*1e6; % micron per second
s = 0.4*w; 

%% Seq

Seq = A1*w/(alpha2*sqrt(rm));


disp(strcat('Critical S:',num2str(100*(s_crit-1))))
disp(strcat('Equilibrium S:',num2str(100*Seq)))
disp(strcat('Critical r2:',num2str(r2_crit)))
disp(strcat('Steady r2:',num2str(rm)))

epsilon = A1*s/(alpha2*sqrt(rm));


%%

time_vals = 10.^(linspace(1,4,20));
Num_ens = 1e2;

%% Initialise solution
dt = 1e-2;
frac_act = zeros(1,length(time_vals));

for K = 1:length(time_vals)

aux_fa = 0;

tmax = time_vals(K);
n = floor(tmax/dt);

for nens = 1:Num_ens    
    R = zeros(1,n);
    S = zeros(1,n);
    R(1) = rm;
    
    for i=1:n-1
        S(i+1) = S(i) - dt*alpha2*S(i)*sqrt(R(i)) + dt*A1*w + sqrt(dt)*A1*s*randn();
        R(i+1) = max(0,R(i) + dt*2*A3*(S(i)+1-RHeq(sqrt(R(i)),k)));
        if R(i+1)>r2_crit
            aux_fa = aux_fa + 1;
            break
        end
    end
    m_ss_1 = m_ss_1 + mean(S(1:i));
end
frac_act(K) = aux_fa/Num_ens;
K/length(time_vals)
end
%%

save('frac_activation_s4_full.mat','frac_act')

