
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
s = 0.3*w; 

%% Seq

Seq = A1*w/(alpha2*sqrt(rm));


disp(strcat('Critical S:',num2str(100*(s_crit-1))))
disp(strcat('Equilibrium S:',num2str(100*Seq)))
disp(strcat('Critical r2:',num2str(r2_crit)))
disp(strcat('Steady r2:',num2str(rm)))

epsilon = A1*s/(alpha2*sqrt(rm));

%%

tmax = 1e3;
dt = 1e-2;
n = floor(tmax/dt);


%% Initialise solution

R = zeros(1,n);
R(1) = rm;

R_eq1 = zeros(1,n);
R_eq2 = zeros(1,n);
options = optimset('TolFun',1e-15,'Display','off');


for i=1:n-1
    
    R(i+1) = max(0,R(i) + dt*2*A3*(Seq+1-RHeq(sqrt(R(i)),k))+sqrt(dt)*A3*epsilon*randn()  );

end


%%


figure
plot(dt*[1:n],R)
hold on
plot(dt*[1:n],r2_crit*ones(1,n),'color','red')

