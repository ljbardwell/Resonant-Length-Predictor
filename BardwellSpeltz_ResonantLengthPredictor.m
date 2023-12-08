% Lydia Bardwell Speltz
% 12-8-23
% Resonant length prediction assuming constant electric field 

%% Inputs to program 
clear all;
clc;
close all;

fprintf('This resonant length prediction tool is from the following reference: Bardwell Speltz L, Lee SK, Shu Y, Tarasek M, Trzasko J, Foo TK, Bernstein MA. Modeling and measurement of lead tip heating and resonant length for implanted, insulated wires.\n');


w=2*pi()*63.87*10^6; % frequency in rad/s (2*pi*f)(1.5T)
%w=2*pi()*127.74*10^6;  % frequency in rad/s (2*pi*f) (3T)
mu=4*pi()*10^-7; % permeability in H/m
a=0.39/1000; % radius of wire in m
b=0.625/1000; % radius of wire+insulation in m
e0=8.85418782*10^-12; % absolute permittivity F/m
cond=0.47; % conductivity in S/m
perm=80*e0; % permittivity  in F/m
E=1; % constant E field in V/m
er=2.3; % relative permittivity of insulation 
gamma=1; % reflection coeeficient for the transmission line model 

%% Wavenumber of Medium 
kt_real=w.*sqrt((perm.*mu)./2).*((sqrt(1+(cond./(perm.*w)).^2)) + 1 ).^(1/2); 
kt_imag=w.*sqrt((perm.*mu)./2).*((sqrt(1+(cond./(perm.*w)).^2)) - 1 ).^(1/2);
kt = complex(kt_real, kt_imag); 

%% Wavenumber of Insulator
ki_real=w*sqrt((e0*er)*mu); 
ki_imag=0; 
ki = complex(ki_real, ki_imag); 

%% King Wavenumber 
H0 = besselh(0,kt*b);
H1 = besselh(1,kt*b);
F=H0/((kt*b)*H1);

k=ki*(1+(F/(log(b/a))))^(1/2); 
kI=-imag(k);
kR=real(k);
k=complex(kR,kI);

%% Simple Exponential Model 
d=linspace(0,50,1001)./100; %Length of wire in m
% Voltage Equation 
V=(E/(1i.*k)).*(((exp(-kI.*d)).*cos(kR.*d)-1)+1i.*(exp(-kI.*d)).*sin(kR.*d));
%Voltage Squared Equation
Vs=((E^2)/(abs(k).^2)).*(1+exp(2.*kI.*d)-2.*exp(kI.*d).*cos(kR.*d));
[Y,I]=max(Vs); 

% Derivative of voltage squared equation
dVs=(kI.*exp(kI.*d)-kI.*cos(kR.*d)+kR.*sin(kR.*d));

figure; plot(d*100,Vs);ylabel('Voltage Squared(V^2)'); xlabel('Length (cm)'); title('SEM');
hold on;
plot(d(1,I)*100,Vs(1,I),'o'); legend('Voltage Squared','Resonant Length');
clear Vs V

% Find where the derivative is equal to 0
myfun = @(d,E,k,kI,kR) (((E^2)/abs(k)^2).*(-2.*kI.*exp(2.*kI.*d)+2.*kI.*exp(kI.*d).*cos(kR.*d)+2.*kR.*exp(kI.*d).*sin(kR.*d)));  % parameterized function
fun = @(z) myfun(z,E,k,kI,kR);    
x = fzero(fun,0.24);

fprintf('The resonant length using the King Wavenumber in the SEM is %0.4f cm\n',d(1,I)*100);

%% Transmission Line Model 
d=linspace(0,50,1001)./100; %Length of wire in m

%Voltage Squared Equation
Vs=E.^2.*((1./(exp(-2.*kI.*d)+(gamma.^2).*exp(2.*kI.*d)-(2.*gamma.*cos(2.*kR.*d))))).*((1./abs(k.^2)).*((gamma.^2).*exp(2.*kI.*d)+(2.*gamma.*cos(2.*kR.*d))+((1+gamma).^2)+(exp(-2.*kI.*d))-(2.*cos(kR.*d).*(1+gamma).*((gamma.*exp(kI.*d))+exp(-kI.*d)))));

% Find where the voltage is a max
[Y,I]=max(Vs); 

figure; plot(d*100,Vs);ylabel('Voltage Squared(V^2)'); xlabel('Length (cm)'); title('TLM');
hold on;
plot(d(1,I)*100,Vs(1,I),'o'); legend('Voltage Squared','Resonant Length');

fprintf('The resonant length using the King Wavenumber in the TLM is %0.4f cm\n',d(1,I)*100);
