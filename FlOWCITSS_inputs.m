%% User Inputs 
tIncrement=0.002;
tMax=70;
linspaceT=0:tIncrement:tMax;
lengthT=length(linspaceT);
fileNameAddedMassRadDamp='rm6_1.txt';
fileNameExcForce='rm6_3.txt';
fileNameZVel='rm6_vz.txt';
fileNameHst='rm6_hst.txt';
freqValIncrement=0.01;
reffreqValIncrement=0.005;
freqValMin=0;
freqValMax=3.5;
betaValIncrement=5;
refbetaValIncrement=0.1;
betaValMax=180;
betaValMin=0;
dof=3;
dofValues=[1 3 5];
ULEN=37.7;
g=9.80665;
rho=1000;
n=[0 0 1];
internalArea=27*17.5;
numField=231;
heading=180; 
xFieldPtStart=-13.86;
yFieldPtStart=-13.49;
xFieldPtInc=1.748;
yFieldPtInc=1.349;
xFieldPtEnd=3.62;
yFieldPtEnd=13.49;

% SIMULINK VARIABLES
hydrostatics=importdata(fileNameHst,' ');
Mm=1/(2.0270*10^6); %1/device mass
I5=1/(4.1625*10^8);
p0=101325; %atm pressure in Pascals
rhoa= 1.225; %air density
gamma=1.4; %spec. heat ratio
V0=27*17.5*10; %average air chamber volume

Ru=20; %load resistance
Re=10;
%1st dof
A1=472.5;
Cd1=2.5;
SbS1=+488000; %1mode hydrostatic and external device stiffness
Rq1=rho*A1*Cd1; %1mode quadratic damping coefficient
Rc1= 0; % 1mode coulomb damping coefficient
r1=internalArea*0; %transformation factor 1st dof
Red1= 91000; %external damping
%3rd dof
A3=945;
Cd3=5;
SbS3=hydrostatics.data(15,3)*(rho*g*(ULEN^2))+11100; %1mode hydrostatic and external device stiffness
SbS35=hydrostatics.data(17,3)*(rho*g*(ULEN^3));
Rq3=rho*A3*Cd3; %1mode quadratic damping coefficient
% Rc3= 0; % 1mode coulomb damping coefficient
r3=internalArea*1; %transformation factor 1st dof
Red3= 619000; %external damping
%5th dof
SbS5=hydrostatics.data(29,3)*(rho*g*(ULEN^4))+9920000; %1mode hydrostatic and external device stiffness
SbS53=hydrostatics.data(27,3)*(rho*g*(ULEN^3));
Rq5=0; %1mode quadratic damping coefficient
% Rc5= 0; % 1mode coulomb damping coefficient
r5=-internalArea*(5); %transformation factor 5th dof, multiplied by x coordinate of internal free surface
Red5=190000000; %external damping


% MOORING VARIABLES
mooringLineLength=810;
mooringLineDepth=400;
mooringLineX1=600;
lineWeight=74*9.81; %N/m of chain due to gravity
iterationThreshold= 0.0025;
startingT0=1000;
alphai0=lineWeight/startingT0;
Tv0=5.81*10^5+148;

% WAVE ENVIRONMENT 
tIncrement=0.001;
linspaceT=0:tIncrement:100;
syms t waveNum
waveHeight=3;
Period=10;
omegaFreq=2*pi/Period;
avgBodyDepth=-8.5;



