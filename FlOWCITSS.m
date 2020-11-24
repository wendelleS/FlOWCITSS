%% FlOWCITSS
%Written by Wendelle Sparrer at Virginia Tech
%November 24th, 2020

 
%% Establish Preliminary Values
load('FlOWCITSS_inputs.m')
freqVal=freqValMin:freqValIncrement:freqValMax;
betaVal=betaValMin:betaValIncrement:betaValMax;
refinedFreqVal=0:0.005:freqValMax;
headingIndex=find(betaVal==heading);
numFreq=length(freqVal);
numBeta=length(betaVal);
xFieldPtVal=xFieldPtStart:xFieldPtInc:xFieldPtEnd;
yFieldPtVal=yFieldPtStart:yFieldPtInc:yFieldPtEnd;
xValLength=length(xFieldPtVal);
yValLength=length(yFieldPtVal);
tempBetaVal=betaValMin:refbetaValIncrement:betaValMax;
%% Import, Sort, and Dimensionalize Z velocity Data from WAMIT 
%sort the z velocity import
vz=importdata(fileNameZVel,' ');
freeSurfVelocityZr=vz.data(1:numField,1:(2+2*dof));
freeSurfVelocityZd=vz.data((numField+1):numField+numField*numBeta,1:5);
vz.data(1:(numField+numField*numBeta),:)=[];
for l=1:(numFreq)
    freeSurfVelocityZr=cat(1,freeSurfVelocityZr,vz.data(1:numField,1:(2+2*dof)));
    freeSurfVelocityZd=cat(1,freeSurfVelocityZd,vz.data((numField+1):(numField+numField*numBeta),1:5));
    vz.data(1:(numField+numField*numBeta),:)=[];
end

%% Evaluate for the Excitation Volume Flow FRF (q)
q=zeros(numBeta*numFreq,3);
for k=1:(numFreq+1)
    for l=1:numBeta
        %establish frequency and headings
        q((k-1)*numBeta+l,1:2)=freeSurfVelocityZd((((k-1)*(numBeta)+l)-1)*numField+1,1:2);
  
        % evaluate real and imag parts 
         realGrid=zeros(yValLength,xValLength);
         for r=1:yValLength
             realGrid(r,:)=transpose(freeSurfVelocityZd((((k-1)*(numBeta)+l)-1)*numField+(r-1)*xValLength+1:(((k-1)*(numBeta)+l)-1)*numField+(r*xValLength),4));
         end
          imagGrid=zeros(yValLength,xValLength);
         for r=1:yValLength
             imagGrid(r,:)=transpose(freeSurfVelocityZd((((k-1)*(numBeta)+l)-1)*numField+(r-1)*xValLength+1:(((k-1)*(numBeta)+l)-1)*numField+(r*xValLength),5));
         end
         realsum=cumtrapz(yFieldPtVal,cumtrapz(xFieldPtVal,realGrid,2));
         realpt=realsum(yValLength,xValLength);
         
         imagsum=cumtrapz(yFieldPtVal,cumtrapz(xFieldPtVal,imagGrid,2));
         imagpt=imagsum(yValLength,xValLength);
         
       % Dimensionalize and add together real and imaginary parts of q
         if q((k-1)*numBeta+l,1)==0||q((k-1)*numBeta+l,1)==-1
            q((k-1)*numBeta+l,3)=realpt+1i*imagpt;
         else
            q((k-1)*numBeta+l,3)=(1i.*g.*(realpt+1i*imagpt))./(ULEN.*q((k-1)*numBeta+l,1));
         end
         
    end
end
%sort in order of wave heading, then frequency

q=sortrows(q);
qInfIndex=find(q(:,1)==-1);
qInf=q(qInfIndex(1):qInfIndex(numBeta),:);
q(qInfIndex(1):qInfIndex(numBeta),:)=[ ];
q=sortrows(q,[2,1]);
% %extract the desired wave heading
qHeading=q((headingIndex*(numFreq)-numFreq+1):headingIndex*(numFreq),:);
% 
% Plot some of the results
figure
plot(qHeading(:,1),abs(qHeading(:,3)),'k')
xlabel('Frequency [rad/s]')
ylabel('m^3/s')
title('Volumetric Flow FRF')

%% Evaluate for the Radiation Admittance FRF (Y) Read and Imaginary Components (G and B) 

syms waveNum
q=sortrows(q,[1,2]);
G=zeros(numFreq,2);

for k=1:numFreq
    %obtain frequency value
    G(k,1)=q((k-1)*numBeta+1,1);
 
    %integrate across beta
    interpolatedq=interp1(q((k-1)*numBeta+1:(k-1)*numBeta+numBeta,2),q((k-1)*numBeta+1:(k-1)*numBeta+numBeta,3), tempBetaVal);
    radians=deg2rad(tempBetaVal);
    G(k,2)=trapz(radians,(abs(interpolatedq).^2));
    %multiply by the appropriate coeffs
    G(k,2)=G(k,2)*(((G(k,1)^3))/(2*pi*rho*g^3));

end
G(1,2)=0;

G=[transpose(refinedFreqVal),transpose(interp1(G(:,1),G(:,2),refinedFreqVal))];

maxFreqInd=length(refinedFreqVal);
B=zeros(maxFreqInd-2,2);
 for k=3:(maxFreqInd-2)
    B(k,1)=G(k,1);
    beforeFreq=trapz(G(1:(k-1),1),(G(1:(k-1),2)./(B(k,1)^2-G(1:(k-1),1).^2)));
    afterFreq=trapz(G((k+1):maxFreqInd,1),(G((k+1):maxFreqInd,2)./(B(k,1)^2-G((k+1):maxFreqInd,1).^2)));
    B(k,2)=beforeFreq+afterFreq;
    B(k,2)=(-2.*B(k,1)/pi).*B(k,2);
 end
 B(1,:)=[];

figure
plot(G(:,1),G(:,2),'k',B(:,1),B(:,2),'-r','LineWidth', 1.5)
legend('G, Re[Y]','B, Im[Y]')
xlabel('Frequency [rad/s]')
ylabel('Radiation Admittance [m^4*s/kg]')


%% Evaluate for the Coupling FRF (H)
H=zeros(numFreq,1+dof);

for k=1:(numFreq+1)
    for j=1:dof
        %establish frequency
        H(k,1)=freeSurfVelocityZr((k*numField-numField+1),1);
        % evaluate H, real and imaginary parts  
        realGrid=zeros(yValLength,xValLength);
         for r=1:length(yFieldPtVal)
             realGrid(r,:)=transpose(freeSurfVelocityZr((k*numField-numField+1+(r-1)*xValLength):(k*numField-numField+r*xValLength),(2+2*j-1)));
         end
          imagGrid=zeros(yValLength,xValLength);
         for r=1:length(yFieldPtVal)
             imagGrid(r,:)=transpose(freeSurfVelocityZr((k*numField-numField+1+(r-1)*xValLength):(k*numField-numField+r*xValLength),(2+2*j)));
         end
         realsum=cumtrapz(yFieldPtVal,cumtrapz(xFieldPtVal,realGrid,2));
         realpt=realsum(yValLength,xValLength);
         imagsum=cumtrapz(yFieldPtVal,cumtrapz(xFieldPtVal,imagGrid,2));
         imagpt=imagsum(yValLength,xValLength);
% non-dimensionalize H
        if H(k,1)==0||H(k,1)==-1
            H(k,1+j)=((ULEN^(n(j)))*(-realpt-1i*imagpt));
        else
            H(k,1+j)=(g.*ULEN^(n(j))*(-realpt-1i*imagpt))/((ULEN)*(H(k,1).^2));
        end
    end
end

% extract out C(inf)
H=sortrows(H);
HInfIndex=find(H(:,1)==-1);
cInf=H(HInfIndex,:);
H(HInfIndex,:)=[ ];
H_J=[H(:,1),imag(H(:,2)),imag(H(:,3)),imag(H(:,4))];
H_C=[H(:,1),real(H(:,2)),real(H(:,3)),real(H(:,4))];
C1inf= cInf(1,2); %mode 1 infinite C
C3inf= cInf(1,3); %mode 3 infinite C
C5inf= cInf(1,4); %mode 5 infinite C

figure
plot(H_J(:,1),H_J(:,3),'k',H_C(:,1),H_C(:,3),':r','LineWidth', 1.5)
legend('J, Im[H]','C, Re[H]')
xlabel('Frequency [rad/s]')
ylabel('[m^2]')
title('Coupling Term Heave FRFs')
% 
figure
plot(H_J(:,1),H_J(:,4),'k',H_C(:,1),H_C(:,4),':r','LineWidth', 1.5)
legend('J, Im[H]','C, Re[H]')
xlabel('Frequency [rad/s]')
ylabel('[m^2]')
title('Coupling Term Pitch FRFs')
%% import the rest of the pertinent WAMIT Files and clear other data 
% clear freeSurfVelocityZr freeSurfVelocityZd
AB=readtable(fileNameAddedMassRadDamp);
excForce=readtable(fileNameExcForce);
%% Evaluate for the Excitation Force FRF (f)
excForce=sortrows(excForce,[3 2 1]);
% Dimensionalize
m=zeros((numFreq-1)*numBeta*dof,2);
for k=1:dof
    if dofValues(k)<4
        m((k-1)*(numFreq-1)*numBeta+1:k*(numFreq-1)*numBeta,1:2)=2*ones((numFreq-1)*numBeta,2);
    else
        m((k-1)*(numFreq-1)*numBeta+1:k*(numFreq-1)*numBeta,1:2)=3*ones((numFreq-1)*numBeta,2);
    end
end
excForce{:,6:7}=excForce{:,6:7}.*rho.*g.*(ULEN.^m);
excForce=sortrows(excForce, [2,3,1]);
excForce1=excForce{1:(numFreq-1),:};
excForce3=excForce{numFreq:(numFreq-1)*2,:};
excForce5=excForce{numFreq+numFreq-1:(numFreq-1)*3,:};


figure
plot(excForce3(:,1),abs(excForce3(:,6)+1i.*excForce3(:,7)),'k','LineWidth', 1.5)
xlabel('Frequency [rad/s]')
ylabel('m^3/s')
title('f3')
figure
plot(excForce5(:,1),abs(excForce5(:,6)+1i.*excForce5(:,7)),'k','LineWidth', 1.5)
xlabel('Frequency [rad/s]')
ylabel('m^3/s')
title('f5')
figure
plot(excForce1(:,1),abs(excForce1(:,6)+1i.*excForce1(:,7)),'k','LineWidth', 1.5)
xlabel('Frequency [rad/s]')
ylabel('m^3/s')
title('f1')

%% Evaluate for the Radiation Impedance FRF (Z)

AB=sortrows(AB,[2,3,1]);

numFreq=numFreq+1;
m=zeros(((numFreq)*dof*dof),1);
for k=1:dof^2
    if (AB{(k-1)*numFreq+1,2}<4) && (AB{(k-1)*numFreq+1,3}<4)
        m((k-1)*numFreq+1:k*numFreq,1)=3*ones(numFreq,1);
    elseif (AB{(k-1)*numFreq+1,2}<4) && (AB{(k-1)*numFreq+1,3}>3)
        m((k-1)*numFreq+1:k*numFreq,1)=4*ones(numFreq,1);
    elseif (AB{(k-1)*numFreq+1,2}>3) && (AB{(k-1)*numFreq+1,3}<4)
        m((k-1)*numFreq+1:k*numFreq,1)=4*ones(numFreq,1);
    elseif (AB{(k-1)*numFreq+1,2}>3) && (AB{(k-1)*numFreq+1,3}>3)
        m((k-1)*numFreq+1:k*numFreq,1)=5*ones(numFreq,1);
    end
end

AB{:,4}=AB{:,4}.*rho.*(ULEN.^m);
AB{:,5}=AB{:,5}.*AB{:,1}.*rho.*(ULEN.^m);
AB=sortrows(AB,[1 2 3]);
addedMassInf=AB{1:dof^2,:};
AB=AB((dof^2+1):dof^2*numFreq,:);
AB=sortrows(AB,[2,3,1]);
numFreq=numFreq-1;

b11=[AB{1:numFreq,1},AB{1:numFreq,5}];
a11=[AB{1:numFreq,1},AB{1:numFreq,4}];
AB=AB(numFreq+1:height(AB),:);
b13=[AB{1:numFreq,1},AB{1:numFreq,5}];
a13=[AB{1:numFreq,1},AB{1:numFreq,4}];
AB=AB(numFreq+1:height(AB),:);
b15=[AB{1:numFreq,1},AB{1:numFreq,5}];
a15=[AB{1:numFreq,1},AB{1:numFreq,4}];
AB=AB(numFreq+1:height(AB),:);

b31=[AB{1:numFreq,1},AB{1:numFreq,5}];
a31=[AB{1:numFreq,1},AB{1:numFreq,4}];
AB=AB(numFreq+1:height(AB),:);
b33=[AB{1:numFreq,1},AB{1:numFreq,5}];
a33=[AB{1:numFreq,1},AB{1:numFreq,4}];
AB=AB(numFreq+1:height(AB),:);
b35=[AB{1:numFreq,1},AB{1:numFreq,5}];
a35=[AB{1:numFreq,1},AB{1:numFreq,4}];
AB=AB(numFreq+1:height(AB),:);


b51=[AB{1:numFreq,1},AB{1:numFreq,5}];
a51=[AB{1:numFreq,1},AB{1:numFreq,4}];
AB=AB(numFreq+1:height(AB),:);
b53=[AB{1:numFreq,1},AB{1:numFreq,5}];
a53=[AB{1:numFreq,1},AB{1:numFreq,4}];
AB=AB(numFreq+1:height(AB),:);
b55=[AB{1:numFreq,1},AB{1:numFreq,5}];
a55=[AB{1:numFreq,1},AB{1:numFreq,4}];
AB=AB(numFreq+1:height(AB),:);


a11Inf=addedMassInf(1,4);
a13Inf=addedMassInf(2,4);
a15Inf=addedMassInf(3,4);
a31Inf=addedMassInf(4,4);
a33Inf=addedMassInf(5,4);
a35Inf=addedMassInf(6,4);
a51Inf=addedMassInf(7,4);
a53Inf=addedMassInf(8,4);
a55Inf=addedMassInf(9,4);


figure
plot(a31(:,1),a31(:,2),'k',a33(:,1),a33(:,2),':r', a35(:,1),a35(:,2),'-.m','LineWidth', 1.5)
legend('31','33', '35')
xlabel('Frequency [rad/s]')
ylabel('added mass [kg]')
title('added mass FRFs')

figure
plot(b11(:,1),b11(:,2),'k',b13(:,1),b13(:,2),':r', b15(:,1),b15(:,2),'-.m','LineWidth', 1.5)
legend('11','13', '15')
xlabel('Frequency [rad/s]')
ylabel('Radiation damping [-]')
title('Radiation damping FRFs')

figure
plot(b31(:,1),b31(:,2),'k',b33(:,1),b33(:,2),':r', b35(:,1),b35(:,2),'-.m','LineWidth', 1.5)
legend('31','33', '35')
xlabel('Frequency [rad/s]')
ylabel('Radiation damping [-]')
title('Radiation damping FRFs')
figure
plot(b51(:,1),b51(:,2),'k',b53(:,1),b53(:,2),':r', b55(:,1),b55(:,2),'-.m','LineWidth', 1.5)
legend('51','53', '55')
xlabel('Frequency [rad/s]')
ylabel('Radiation damping [-]')
title('Radiation damping FRFs')
%% State Space Realization
opts.cmplx_ss=0;


%% 11 Rad Impedance
%First, generate the proper FRF for realization
b11(1,2)=0;
K11(:,1)=b11(:,1);
K11(:,2)=b11(:,2)+1i.*b11(:,1).*(a11(:,2)-a11(numFreq,2));
%Next, apply vecfit3
Ns=length(K11);
f=transpose(K11(:,2));
s=transpose(1i.*K11(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=10;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AK11=full(SER.A);
BK11=SER.B;
CK11=SER.C;
DK11=SER.D;
EK11=SER.E;
K11_vecfit=ss(AK11,BK11,CK11,DK11);


%% 13 Rad Impedance
%First, generate the proper FRF for realization
b13(1,2)=0;
K13(:,1)=b13(:,1);
K13(:,2)=b13(:,2)+1i.*b13(:,1).*(a13(:,2)-a13(numFreq,2));
%Next, apply vecfit3
Ns=length(K13);
f=transpose(K13(:,2));
s=transpose(1i.*K13(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=10;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AK13=full(SER.A);
BK13=SER.B;
CK13=SER.C;
DK13=SER.D;
EK13=SER.E;
K13_vecfit=ss(AK13,BK13,CK13,DK13);


%% 15 Rad Impedance
%First, generate the proper FRF for realization
b15(1,2)=0;
K15(:,1)=b15(:,1);
K15(:,2)=b15(:,2)+1i.*b15(:,1).*(a15(:,2)-a15(numFreq,2));
%Next, apply vecfit3
Ns=length(K15);
f=transpose(K15(:,2));
s=transpose(1i.*K15(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=10;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AK15=full(SER.A);
BK15=SER.B;
CK15=SER.C;
DK15=SER.D;
EK15=SER.E;
K15_vecfit=ss(AK15,BK15,CK15,DK15);

%% 31 Rad Impedance
%First, generate the proper FRF for realization
b31(1,2)=0;
K31(:,1)=b31(:,1);
K31(:,2)=b31(:,2)+1i.*b31(:,1).*(a31(:,2)-a31(numFreq,2));
%Next, apply vecfit3
Ns=length(K31);
f=transpose(K31(:,2));
s=transpose(1i.*K31(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=10;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AK31=full(SER.A);
BK31=SER.B;
CK31=SER.C;
DK31=SER.D;
EK31=SER.E;
K31_vecfit=ss(AK31,BK31,CK31,DK31);


%% 33 Rad Impedance
%First, generate the proper FRF for realization
b33(1,2)=0;
K33(:,1)=b33(:,1);
K33(:,2)=b33(:,2)+1i.*b33(:,1).*(a33(:,2)-a33(numFreq,2));
%Next, apply vecfit3
Ns=length(K33);
f=transpose(K33(:,2));
s=transpose(1i.*K33(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=10;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AK33=full(SER.A);
BK33=SER.B;
CK33=SER.C;
DK33=SER.D;
EK33=SER.E;
K33_vecfit=ss(AK33,BK33,CK33,DK33);


%% 35 Rad Impedance
%First, generate the proper FRF for realization
b35(1,2)=0;
K35(:,1)=b35(:,1);
K35(:,2)=b35(:,2)+1i.*b35(:,1).*(a35(:,2)-a35(numFreq,2));
%Next, apply vecfit3
Ns=length(K35);
f=transpose(K35(:,2));
s=transpose(1i.*K35(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=10;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AK35=full(SER.A);
BK35=SER.B;
CK35=SER.C;
DK35=SER.D;
EK35=SER.E;
K35_vecfit=ss(AK35,BK35,CK35,DK35);

%% K51
b51(1,2)=0;
K51(:,1)=b51(:,1);
K51(:,2)=b51(:,2)+1i.*b51(:,1).*(a51(:,2)-a51(numFreq,2));
%Next, apply vecfit3
Ns=length(K51);
f=transpose(K51(:,2));
s=transpose(1i.*K51(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=10;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AK51=full(SER.A);
BK51=SER.B;
CK51=SER.C;
DK51=SER.D;
EK51=SER.E;
K51_vecfit=ss(AK51,BK51,CK51,DK51);


%% 53 Rad Impedance
%First, generate the proper FRF for realization
b53(1,2)=0;
K53(:,1)=b53(:,1);
K53(:,2)=b53(:,2)+1i.*b53(:,1).*(a53(:,2)-a53(numFreq,2));
%Next, apply vecfit3
Ns=length(K53);
f=transpose(K53(:,2));
s=transpose(1i.*K53(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=10;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AK53=full(SER.A);
BK53=SER.B;
CK53=SER.C;
DK53=SER.D;
EK53=SER.E;
K53_vecfit=ss(AK53,BK53,CK53,DK53);

%% 55 Rad Impedance
%First, generate the proper FRF for realization
b55(1,2)=0;
K55(:,1)=b55(:,1);
K55(:,2)=b55(:,2)+1i.*b55(:,1).*(a55(:,2)-a53(numFreq,2));
%Next, apply vecfit3
Ns=length(K55);
f=transpose(K55(:,2));
s=transpose(1i.*K55(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=10;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AK55=full(SER.A);
BK55=SER.B;
CK55=SER.C;
DK55=SER.D;
EK55=SER.E;
K55_vecfit=ss(AK55,BK55,CK55,DK55);

%% Radiation Admittance 
%First, generate the proper FRF for realization
clear Y Gss Bss
Gss(:,1)=transpose([0:reffreqValIncrement/2:freqValMax-freqValIncrement]);
Gss(:,2)=interp1(G(:,1),G(:,2),Gss(:,1));
Bss(:,2)=interp1(B(:,1),B(:,2),Gss(:,1));
Y(:,1)=Gss(:,1);
Y(:,2)=Gss(:,2)+1i.*Bss(:,2);

%Next, apply vecfit3
Ns=length(Y);
f=transpose(Y(:,2));
s=transpose(1i.*Y(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax-freqValIncrement,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=12;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AY=full(SER.A);
BY=SER.B;
CY=SER.C;
DY=SER.D;
EY=SER.E;
Y_vecfit=ss(AY,BY,CY,DY);
%% Radiation Coupling 1
% h1Freqs=0:0.01:39.9;
H1(:,1)=H_C(:,1);
H1(:,2)=(H_C(:,2)-H_C(numFreq,2))+1i.*H_J(:,2);
%Next, apply vecfit3
Ns=length(H1);
f=transpose(H1(:,2));
s=transpose(1i.*H1(:,1));
weight=ones(1,Ns);
N=50; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=15;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AH1=full(SER.A);
BH1=SER.B;
CH1=SER.C;
DH1=SER.D;
EH1=SER.E;
H1_vecfit=ss(AH1,BH1,CH1,DH1);
%% Radiation Coupling 3
H3(:,1)=H_C(:,1);
H3(:,2)=(H_C(:,3)-H_C(numFreq,3))+1i.*H_J(:,3);
%Next, apply vecfit3
Ns=length(H3);
f=transpose(H3(:,2));
s=transpose(1i.*H3(:,1));
weight=ones(1,Ns);
N=160; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=15;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AH3=full(SER.A);
BH3=SER.B;
CH3=SER.C;
DH3=SER.D;
EH3=SER.E;
H3_vecfit=ss(AH3,BH3,CH3,DH3);
%% Radiation Coupling 5
H5(:,1)=H_C(:,1);
H5(:,2)=(H_C(:,4)-H_C(numFreq,4))+1i.*H_J(:,4);
%Next, apply vecfit3
Ns=length(H5);
f=transpose(H5(:,2));
s=transpose(1i.*H5(:,1));
weight=ones(1,Ns);
N=120; %Order of approximation 
%Complex conjugate pairs, linearly spaced: 
bet=linspace(0,freqValMax,N/2); 
poles=[]; 
for n=1:length(bet)   
    alf=-bet(n)*1e-2;   
    poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
end   
Niter=10;
for iter =1:Niter    
    [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);
end 
AH5=full(SER.A);
BH5=SER.B;
CH5=SER.C;
DH5=SER.D;
EH5=SER.E;
H5_vecfit=ss(AH5,BH5,CH5,DH5);


%% Excitation Inputs from Test Data 

excForce1=[0,0,0,0,0,excForce1(1,6),excForce1(1,7);excForce1];
excForce1r=transpose(interp1(excForce1(:,1),excForce1(:,6),refinedFreqVal));
excForce1i=transpose(interp1(excForce1(:,1),excForce1(:,7),refinedFreqVal));
excForce1=[transpose(refinedFreqVal),excForce1r,excForce1i];

excForce3=[0,0,0,0,0,excForce3(1,6),excForce3(1,7);excForce3];
excForce3r=interp1(excForce3(:,1),excForce3(:,6),refinedFreqVal);
excForce3i=interp1(excForce3(:,1),excForce3(:,7),refinedFreqVal);
excForce3=[transpose(refinedFreqVal),transpose(excForce3r),transpose(excForce3i)];

excForce5=[0,0,0,0,0,excForce5(1,6),excForce5(1,7);excForce5];
excForce5r=interp1(excForce5(:,1),excForce5(:,6),refinedFreqVal);
excForce5i=interp1(excForce5(:,1),excForce5(:,7),refinedFreqVal);
excForce5=[transpose(refinedFreqVal),transpose(excForce5r),transpose(excForce5i)];

excQ=interp1(qHeading(:,1),qHeading(:,3),refinedFreqVal);
excQ=[transpose(refinedFreqVal),transpose(excQ)];

%% create theoretical inputs 
waveAmp=waveHeight/2;
waveNumber=(omegaFreq^2)/g;
propagationDistance=0;
omegaFreq=2*pi/Period;
blendIn=@(t) (1/5)*t;
waveElevation= piecewise(t<5, (1/5)*t*real(waveAmp*exp(1i*omegaFreq*t)), t>=5, real(waveAmp*exp(1i*omegaFreq*t)));
figure
fplot(waveElevation, [0 100])
title('theoretical wave')

blendInVals=double(blendIn(0:tIncrement:5));
blendInValues=[blendInVals,ones(1,length((5+tIncrement):tIncrement:100))];
verticalParticleVelocity=@(t) (pi*waveHeight/Period)*exp(waveNumber*avgBodyDepth)*sin(waveNumber*propagationDistance-omegaFreq*t);
testvVel=double(verticalParticleVelocity(linspaceT));
testvVel=blendInValues.*testvVel;
verticalWaveVelocity=timeseries(transpose(testvVel),transpose(linspaceT));

hParticleVelocity=@(t)  (pi*waveHeight/Period)*exp(waveNumber*avgBodyDepth)*cos(waveNumber*propagationDistance-omegaFreq*t);
testhVel=double(-hParticleVelocity(linspaceT));
testhVel=blendInValues.*testhVel;
horizontalWaveVelocity=timeseries(transpose(testhVel),transpose(linspaceT));

% Create excitation inputs from theoretical values 
f1Real=interp1(excForce1(:,1),excForce1(:,2), omegaFreq);
f1Imag=interp1(excForce1(:,1),excForce1(:,3), omegaFreq);
f1ExcT=@(t) waveAmp*real((f1Real+1i*f1Imag)*exp(1i*omegaFreq*t));
f1Samples=double(f1ExcT(linspaceT));
f1Samples=blendInValues.*f1Samples;
f1Exc=timeseries(transpose(f1Samples),transpose(linspaceT));

f3Real=interp1(excForce3(:,1),excForce3(:,2), omegaFreq);
f3Imag=interp1(excForce3(:,1),excForce3(:,3), omegaFreq);
f3ExcT=@(t) waveAmp*real((f3Real+1i*f3Imag)*exp(1i*omegaFreq*t));
f3Samples=double(f3ExcT(linspaceT));
f3Samples=blendInValues.*f3Samples;
f3Exc=timeseries(transpose(f3Samples),transpose(linspaceT));


f5Real=interp1(excForce5(:,1),excForce5(:,2), omegaFreq);
f5Imag=interp1(excForce5(:,1),excForce5(:,3), omegaFreq);
f5ExcT=@(t) waveAmp*real((f5Real+1i*f5Imag)*exp(1i*omegaFreq*t));
f5Samples=double(f5ExcT(linspaceT));
f5Samples=blendInValues.*f5Samples;
f5Exc=timeseries(transpose(f5Samples),transpose(linspaceT));

qInterpolated=interp1(excQ(:,1),excQ(:,2), omegaFreq);
qImag=imag(qInterpolated);
qReal=real(qInterpolated);
qExcT=@(t) waveAmp*real((qReal+1i*qImag)*exp(1i*omegaFreq*t));
qSamples=double(qExcT(linspaceT));
qSamples=blendInValues.*qSamples;
qExc=timeseries(transpose(qSamples),transpose(linspaceT));

