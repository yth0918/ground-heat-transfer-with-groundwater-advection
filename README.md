function fx=double %the mainprogram of the coupled model for the outlet water temperature 
clear;  %
clc;  %
options=optimset('display','off');
tx=0:3600:3600*24;  % the time and the time interval
N=length(tx);  %  the number of time vectors, and can be confirmed the solve by "for cyclic" sentence.
for kk=1:1:N
t=tx(kk);  %  assign the value of "KK" to "t", such that "t" is a numerical  
fx(kk)=fminsearch(@(k)zuankong(k,t),[30],options); % Optimization program with parameters
end
%%  
function f=zuankong(k,t) % subroutine of the  model inside borehole of double-U 
G=0;
Tin=35; %the inlet water temperature of BHE
%% Related parameters in double U-tube BHE
T0=17; %initial soil temperature 
H=100; %the borehole length
dz=0.11; % the borehole diameter
LAMTB=1.86;   %  thermal conductivity of grout
do=0.025;      % U-tube outer diameter
di=0.020;      % U-tube inner diameter
D=0.035;          % spacing between two legs of U-tube
LAMTP=0.45;%  thermal conductivity of tube
cp=4200;% Specific heat capacity of water in tube
p=1000;% the density of water
M=0.46;%the mass flow rate of water
K=3547; % Convection coefficient of the water：hi=0.023*Re^0.8*Pr^n*lambdaf/(2*ri);
%%
%%  the effective thermal conductivity of the soils
% lambdas=0.98;%sand gravel
lambdas=3.56;%limestone
% lambdas=1.02;%coarse sand
% lambdas=1.03;%fine sand
% lambdas=2.07;%silt
% lambdas=4.5;%sandstone
%%
%%the expression equation of thermal resistance of double U-tube 
Rp=(log(do/di))/(2*pi*LAMTP)+1/(pi*di*K);       % thermal resistance between the water and tube wall
R11=((log(dz/do)-(((LAMTB-lambdas)*log((dz^2-D^2)/dz^2))/(LAMTB+lambdas)))/(2*pi*LAMTB))+Rp;      
R12=(log(dz/(sqrt(2)*D))-((LAMTB-lambdas)*log((dz^4+D^4)/dz^4))/(2*(LAMTB+lambdas)))/(2*pi*LAMTB);   
R13=(log(dz/(2*D))-((LAMTB-lambdas)*log((dz^2+D^2)/dz^2))/(LAMTB+lambdas))/(2*pi*LAMTB);        
Ra=R11+R13+(2*R12);  %R1'
Rb=(R11^2+R13^2+(2*R11*R13)-(4*R12^2))/R12;  %R12'
Rc=((R11-R13)*(R11^2+R13^2+(2*R11*R13)-(4*R12^2)))/(R13^2+(R11*R13)-(2*R12^2));%R13'

Ra1=(M*cp*Ra)/H;      %R1*
Rb2=(M*cp*Rb)/H;     %R12*
Rc3=(M*cp*Rc)/H;    %R13*
S12=(Rb2*Rc3)/(Rb2+Rc3); % S12
     BAT=sqrt((1/Ra1^2)+(2/(Ra1*S12)));   %S1=R1*     
CHBAT=(exp(BAT)+exp(-BAT))/2;          
SHBAT=(exp(BAT)-exp(-BAT))/2;


%The outside model solves the value of Tb
x=0.055;    % x=the borehole radius 
z=50;  % z is the value of middle section of the borehole  
y=0;
q=M*(Tin-k)*cp/H;  % k is the assumed outlet temperature 

%%the physical property parameter of the soils
%  the physical property parameter of sand gravel
% m(1)=0.98; % the heat conductivity coefficient of sand gravel
% m(2)=0.7e-6;% the thermal diffusivity of sand gravel
% u=3e-5; %the groundwater advection velocity in the sand gracel

%  the physical property parameter of limestone
m(1)=3.56;% the heat conductivity coefficient
m(2)=2.657e-6;% the thermal diffusivity
u=1e-6; %the groundwater advection velocity
%
% % the physical property parameter of coarse sand
% m(1)=1.02; % the heat conductivity coefficient
% m(2)=0.729e-6;% the thermal diffusivity
% u=7.3e-7; %the groundwater advection velocity
% 
% % the physical property parameter of fine sand
% m(1)=1.03;% the heat conductivity coefficient
% m(2)=0.736e-6;% the thermal diffusivity
% u=6.3e-8; %the groundwater advection velocity
% 
% % the physical property parameter of silt
% m(1)=2.07;% the heat conductivity coefficient
% m(2)=0.726e-6;% the thermal diffusivity
% u=1.4e-9; %the groundwater advection velocity
% 
% % the physical property parameter of sandstone
% m(1)=4.5;% the heat conductivity coefficient
% m(2)=1.264e-6;% the thermal diffusivity
% u=4.2e-10; %the groundwater advection velocity

%% subroutine of the  model outside borehole of double-U 
U1=u*(m(1)/m(2))/p/cp;
ds1=1/100;
s1=0:ds1:100;%the interval is 0 to 100(the borehole depth)
F1=0.25./sqrt(x^2+y^2+(z-s1).^2).*(exp(-U1*sqrt(x^2+y^2+(z-s1).^2)./2/m(2)).*erfc((sqrt(x^2+y^2+(z-s1).^2)-U1*t)./2/(sqrt(m(2)*t)))+exp(U1*sqrt(x^2+y^2+(z-s1).^2)./2/m(2)).*erfc((sqrt(x^2+y^2+(z-s1).^2)+U1*t)./2/(sqrt(m(2)*t))));
S1=ds1*cumtrapz(F1);
S1(end);
% 
dl1=1/100;
l1=-100:dl1:0;
E1=0.25./sqrt(x^2+y^2+(z-l1).^2).*(exp(-U1*sqrt(x^2+y^2+(z-l1).^2)./2/m(2)).*erfc((sqrt(x^2+y^2+(z-l1).^2)-U1*t)./2/(sqrt(m(2)*t)))+exp(U1*sqrt(x^2+y^2+(z-l1).^2)./2/m(2)).*erfc((sqrt(x^2+y^2+(z-l1).^2)+U1*t)./2/(sqrt(m(2)*t))));
L1=dl1*cumtrapz(E1);

%% 
L1(end);
T1=(S1(end)-L1(end))*q/2/pi/m(1)*exp(U1*x/2/m(2));

Tb=T1+T0; %the value of borehole wall "Tb" by using the model outside borehole with groundwater advection

Tout=Tb+(Tin-Tb)*((BAT*Ra1*CHBAT)-SHBAT)/((BAT*Ra1*CHBAT)+SHBAT); %the value of the outlet water temperature by using the model inside borehole

G=(Tout-k)^2+G; % Optimal solution formula

f=G;
  
