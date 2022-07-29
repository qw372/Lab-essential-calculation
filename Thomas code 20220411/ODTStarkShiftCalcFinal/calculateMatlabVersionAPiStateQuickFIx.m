clear all;

calculateUncertainties = 1;

%molecule details
laserEnergy = 9398.5;
APiOneHalfEnergyRelToZero = 15072.09;
APiOneHalfEnergy = 0;
BSigmaEnergy = 17267.4465-APiOneHalfEnergyRelToZero;
ADeltaEnergy = 19108-APiOneHalfEnergyRelToZero;
CPiEnergy = 27384.67-APiOneHalfEnergyRelToZero;
DSigmaEnergy = 27773.8-APiOneHalfEnergyRelToZero;
XSigmaEnergy = 0-APiOneHalfEnergyRelToZero;
% FSigmaEnergy = 32823.5;
% GPiEnergy = 34808.9275;
% HEnergy = 42345;
detuningBSigma = laserEnergy-BSigmaEnergy;
detuningCPi = laserEnergy-CPiEnergy;
detuningDSigma = laserEnergy-DSigmaEnergy;
detuningADelta = laserEnergy-ADeltaEnergy;
% detuningFSigma = laserEnergy-FSigmaEnergy;
% detuningGPi = laserEnergy-GPiEnergy;
% detuningH = laserEnergy-HEnergy;
detuningXSigma = -laserEnergy-XSigmaEnergy;%laser energy has a minus sign here since the X-A interaction is strongest when photon being added to field

% detuningAPiOneHalfCoRotate = laserEnergy+APiOneHalfEnergy;
% detuningAPiThreeHalvesCoRotate = laserEnergy+APiThreeHalvesEnergy;
% detuningBSigmaCoRotate = laserEnergy+BSigmaEnergy;
% detuningCPiOneHalfCoRotate = laserEnergy+CPiEnergyOneHalf;
% detuningCPiThreeHalvesCoRotate = laserEnergy+CPiEnergyThreeHalves;
% detuningDSigmaCoRotate = laserEnergy+DSigmaEnergy;
% detuningFSigmaCoRotate = laserEnergy+FSigmaEnergy;
% detuningGPiCoRotate = laserEnergy+GPiEnergy;
% detuningHCoRotate = laserEnergy+HEnergy;

addCoRotate=1;

if addCoRotate==1
    APiOneHalfEnergy = 15075.6122;
    BSigmaEnergy = 17267.4465;
    ADeltaEnergy = 19108;
    CPiEnergy = 27384.67;
    DSigmaEnergy = 27773.8;
    XSigmaEnergy = 0;
end

hbar=1.0546e-34;
c=2.9979e8;
eps0=8.8542e-12;
mksToDebye = 3.3356e-30;
invCmToHz=2.997e10;
debyeFromLifetime = @(lifetime,energy) sqrt(3*hbar/4/lifetime*4*pi*eps0*c^3/(energy*invCmToHz*2*pi)^3)/mksToDebye;
% debyeAPiOneHalf = 6.1923;%exp
% debyeAPiThreeHalves = 6.2215;%exp
% debyeBStates = 4.8709;%exp
debyeXSigma = debyeFromLifetime(24.1e-9-0e-9,-XSigmaEnergy);%exp
if addCoRotate==1
    debyeXSigma = debyeFromLifetime(24.1e-9-0e-9,APiOneHalfEnergy);%exp
end
debyeBStates = 0.210/.3934;
debyeCPi = 1.1/.3934;%de Melo+Ornellas
debyeDSigma = 2.1/.3934;%de Melo+Ornellas
debyeADelta = 2.711/.3934;

% debyeCPi = debyeFromLifetime(66.3e-9,CPiEnergy)*0;%de Melo+Ornellas
% debyeDSigma = debyeFromLifetime(219e-9,DSigmaEnergy)*0;%de Melo+Ornellas
% debyeFSigma = 0.4/.393*0;%de Melo+Ornellas
% debyeGPi = 0.5/.393*0;%de Melo + Ornellas

% debyeAPiOneHalf = 6.2279;%varun
% debyeAPiThreeHalves = 6.2279;%varun
% debyeBStates = 4.9315;%varun

%
% debyeAPiOneHalf = 6.6816;%comp
% debyeAPiThreeHalves = 6.6807;%comp
% debyeBStates = 5.197;%comp

%experimental details
intensity = 2*51/(pi*(40e-6)^2);
%intensity = 2*51.5/(pi*(39.75e-6)^2);
p=1;%ODT polarization

%fundamental constants
e=1.602e-19;
c=2.9979e8;
hbar=1.0546e-34;
eps0=8.8542e-12;
kB=1.3806e-23;
bohrRadius = 5.2918e-11;
cmMinusOneToFrequency = 2.99793e10*2*pi;
debyeFactor = 0.3934;%convert debye to ea0

prefactorTwo = e^2/4*intensity/c/eps0/hbar/cmMinusOneToFrequency*bohrRadius^2*debyeFactor^2;

%A State values, in terms of |J,F,mF>
% XValues = {[1/2,1,-1],[1/2,1,0],[1/2,1,1],[1/2,0,0],[3/2,1,-1],[3/2,1,0],[3/2,1,1],...
%     [3/2,2,-2],[3/2,2,-1],[3/2,2,0],[3/2,2,1],[3/2,2,2]};

AValues = {[1/2,1,-1],[1/2,1,0],[1/2,1,1],[1/2,0,0]};

%%%FIRST do shifts due to X state.  This will involve calculations
%%%over J=1/2, J=3/2

%J=1/2

XValues = {[1/2,1,-1],[1/2,1,0],[1/2,1,1],[1/2,0,0],[3/2,1,-1],[3/2,1,0],[3/2,1,1],...
    [3/2,2,-2],[3/2,2,-1],[3/2,2,0],[3/2,2,1],[3/2,2,2]};

%X state has
%(1/\sqrt(2))[-(1)^(j+1/2)|sigma=1/2,omega=1/2>+|sigma=-1/2,omega=-1/2>]
%APi_1/2 state has (1/\sqrt(2))[|sigma'=-1/2,omega'=1/2>+(-1)^(j'-1/2)|sigma'=-1/2,omega'=-1/2>]

AShiftDueToXSigma= zeros(1,4);

for i=1:length(AShiftDueToXSigma)
    Fp = AValues{i}(2);
    Jp = AValues{i}(1);
    for j=1:length(XValues)
        
        currAngularTerm = calculateAngularTermsAPi(AValues{i},XValues{j},0,1/2,p);
        AShiftDueToXSigma(i)=AShiftDueToXSigma(i)+currAngularTerm^2;
        
    end
end

%%%NEXT do shifts due interactions between APi and other Pi states (C Pi
%%%specifically)

CValues = {[1/2,1,-1],[1/2,1,0],[1/2,1,1],[1/2,0,0],...
    [3/2,1,-1],[3/2,1,0],[3/2,1,1],[3/2,2,-2],[3/2,2,-1],[3/2,2,0],[3/2,2,1],[3/2,2,2]};

%X state has
%(1/\sqrt(2))[-(1)^(j+1/2)|sigma=1/2,omega=1/2>+|sigma=-1/2,omega=-1/2>]
%APi_1/2 state has (1/\sqrt(2))[|sigma'=-1/2,omega'=1/2>+(-1)^(j'-1/2)|sigma'=-1/2,omega'=-1/2>]

AShiftDueToCPiValues = zeros(1,4);

for i=1:length(AShiftDueToCPiValues)
    Fp = AValues{i}(2);
    Jp = AValues{i}(1);
    for j=1:length(CValues)
        currAngularTerm = calculateAngularTermsAPi(AValues{i},CValues{j},1,1/2,p);
        AShiftDueToCPiValues(i)=AShiftDueToCPiValues(i)+currAngularTerm^2;
        
        
    end
end

%%%NEXT do shifts due to A\Delta state.

%N=0 (B sigma ordering: [J,F,mF,N]

ADeltaValues = {[3/2,1,-1],[3/2,1,0],[3/2,1,1],[3/2,2,-2],[3/2,2,-1],[3/2,2,0],[3/2,2,1],[3/2,2,2]};
%X state has
%(1/\sqrt(2))[-(1)^(j+1/2)|sigma=1/2,omega=1/2>+|sigma=-1/2,omega=-1/2>]
%APi_1/2 state has (1/\sqrt(2))[|sigma'=-1/2,omega'=1/2>+(-1)^(j'-1/2)|sigma'=-1/2,omega'=-1/2>]

AShiftDueToADelta = zeros(1,4);

for i=1:length(AShiftDueToADelta)
    Fp = AValues{i}(2);
    Jp = AValues{i}(1);
    for j=1:length(ADeltaValues)
        
        currAngularTerm = calculateAngularTermsAPi(AValues{i},ADeltaValues{j},2,3/2,p);
        AShiftDueToADelta(i)=AShiftDueToADelta(i)+currAngularTerm^2;
    end
end

% totalEnergyShiftStateFactor = (debyeXSigma^2/detuningXSigma).*(AShiftDueToXSigma)+...
%     (debyeBStates^2/detuningBSigma).*(AShiftDueToXSigma)+...
%     (debyeCPi^2/detuningCPi).*(AShiftDueToCPiValues)+...
%     (debyeDSigma^2/detuningDSigma).*(AShiftDueToXSigma)+...
%     (debyeADelta^2/detuningADelta).*(AShiftDueToADelta);

totalEnergyShiftStateFactor = -(debyeXSigma^2.*(1./(-laserEnergy-APiOneHalfEnergy+XSigmaEnergy))).*(AShiftDueToXSigma)-...
    (debyeBStates^2.*(1./(-laserEnergy-APiOneHalfEnergy+BSigmaEnergy))).*(AShiftDueToXSigma)-...
    (debyeCPi^2.*(1./(-laserEnergy-APiOneHalfEnergy+CPiEnergy))).*(AShiftDueToCPiValues)-...
    (debyeDSigma^2.*(1./(-laserEnergy-APiOneHalfEnergy+DSigmaEnergy))).*(AShiftDueToXSigma)-...
    (debyeADelta^2.*(1./(-laserEnergy-APiOneHalfEnergy+ADeltaEnergy))).*(AShiftDueToADelta);

if addCoRotate==1%NOTE: minus signs are correct here, see ODT grimm paper...really the counter-rot detuning should be pos and the signs above should also be negative, see Grimm ODT review
    p=-p;
    AValues = {[1/2,1,-1],[1/2,1,0],[1/2,1,1],[1/2,0,0]};
    
    %%%FIRST do shifts due to X state.  This will involve calculations
    %%%over J=1/2, J=3/2
    
    %J=1/2
    
    XValues = {[1/2,1,-1],[1/2,1,0],[1/2,1,1],[1/2,0,0],[3/2,1,-1],[3/2,1,0],[3/2,1,1],...
        [3/2,2,-2],[3/2,2,-1],[3/2,2,0],[3/2,2,1],[3/2,2,2]};
    
    %X state has
    %(1/\sqrt(2))[-(1)^(j+1/2)|sigma=1/2,omega=1/2>+|sigma=-1/2,omega=-1/2>]
    %APi_1/2 state has (1/\sqrt(2))[|sigma'=-1/2,omega'=1/2>+(-1)^(j'-1/2)|sigma'=-1/2,omega'=-1/2>]
    
    AShiftDueToXSigma= zeros(1,4);
    
    for i=1:length(AShiftDueToXSigma)
        Fp = AValues{i}(2);
        Jp = AValues{i}(1);
        for j=1:length(XValues)
            
            currAngularTerm = calculateAngularTermsAPi(AValues{i},XValues{j},0,1/2,p);
            AShiftDueToXSigma(i)=AShiftDueToXSigma(i)+currAngularTerm^2;
            
        end
    end
    
    %%%NEXT do shifts due interactions between APi and other Pi states (C Pi
    %%%specifically)
    
    CValues = {[1/2,1,-1],[1/2,1,0],[1/2,1,1],[1/2,0,0],...
        [3/2,1,-1],[3/2,1,0],[3/2,1,1],[3/2,2,-2],[3/2,2,-1],[3/2,2,0],[3/2,2,1],[3/2,2,2]};
    
    %X state has
    %(1/\sqrt(2))[-(1)^(j+1/2)|sigma=1/2,omega=1/2>+|sigma=-1/2,omega=-1/2>]
    %APi_1/2 state has (1/\sqrt(2))[|sigma'=-1/2,omega'=1/2>+(-1)^(j'-1/2)|sigma'=-1/2,omega'=-1/2>]
    
    AShiftDueToCPiValues = zeros(1,4);
    
    for i=1:length(AShiftDueToCPiValues)
        Fp = AValues{i}(2);
        Jp = AValues{i}(1);
        for j=1:length(CValues)
            currAngularTerm = calculateAngularTermsAPi(AValues{i},CValues{j},1,1/2,p);
            AShiftDueToCPiValues(i)=AShiftDueToCPiValues(i)+currAngularTerm^2;
            
            
        end
    end
    
    %%%NEXT do shifts due to A\Delta state.
    
    %N=0 (B sigma ordering: [J,F,mF,N]
    
    ADeltaValues = {[3/2,1,-1],[3/2,1,0],[3/2,1,1],[3/2,2,-2],[3/2,2,-1],[3/2,2,0],[3/2,2,1],[3/2,2,2]};
    %X state has
    %(1/\sqrt(2))[-(1)^(j+1/2)|sigma=1/2,omega=1/2>+|sigma=-1/2,omega=-1/2>]
    %APi_1/2 state has (1/\sqrt(2))[|sigma'=-1/2,omega'=1/2>+(-1)^(j'-1/2)|sigma'=-1/2,omega'=-1/2>]
    
    AShiftDueToADelta = zeros(1,4);
    
    for i=1:length(AShiftDueToADelta)
        Fp = AValues{i}(2);
        Jp = AValues{i}(1);
        for j=1:length(ADeltaValues)
            
            currAngularTerm = calculateAngularTermsAPi(AValues{i},ADeltaValues{j},2,3/2,p);
            AShiftDueToADelta(i)=AShiftDueToADelta(i)+currAngularTerm^2;
        end
    end
    
    
    totalEnergyShiftStateFactorCoRotate = -(debyeXSigma^2.*(1./(laserEnergy-APiOneHalfEnergy+XSigmaEnergy))).*(AShiftDueToXSigma)-...
        (debyeBStates^2.*(1./(laserEnergy-APiOneHalfEnergy+BSigmaEnergy))).*(AShiftDueToXSigma)-...
        (debyeCPi^2.*(1./(laserEnergy-APiOneHalfEnergy+CPiEnergy))).*(AShiftDueToCPiValues)-...
        (debyeDSigma^2.*(1./(laserEnergy-APiOneHalfEnergy+DSigmaEnergy))).*(AShiftDueToXSigma)-...
        (debyeADelta^2.*(1./(laserEnergy-APiOneHalfEnergy+ADeltaEnergy))).*(AShiftDueToADelta);
    totalEnergyShiftStateFactor = totalEnergyShiftStateFactor+totalEnergyShiftStateFactorCoRotate;
end

totalEnergyShiftInMilliKelvin2 = prefactorTwo*totalEnergyShiftStateFactor/kB/1e-3;
totalEnergyShiftInMHz = kB.*totalEnergyShiftInMilliKelvin2.*1e-3./hbar./2./pi./1e6;
alphaT = (totalEnergyShiftStateFactor(4)-totalEnergyShiftStateFactor(2))/(-totalEnergyShiftStateFactor(4));
alphaV = (totalEnergyShiftStateFactor(3)-totalEnergyShiftStateFactor(1))/(2*totalEnergyShiftStateFactor(4));
% if calculateUncertainties ==1
%     aPiOneHalfLifetimeUnc = -2/24.1*1;
%     
%     energiesShiftUpOneHalfLifetimeUnc =  (debyeXSigma^2/detuningXSigma).*(AShiftDueToXSigma)*(1+aPiOneHalfLifetimeUnc);
%     normalEnergyOneHalf = (debyeXSigma^2/detuningXSigma).*(AShiftDueToXSigma)*(1);
%     uncertaintyInShiftRelToScalarOneHalf = (normalEnergyOneHalf-normalEnergyOneHalf(4))-(energiesShiftUpOneHalfLifetimeUnc-energiesShiftUpOneHalfLifetimeUnc(4));
%     uncertaintyInShiftTotalOneHalf = (normalEnergyOneHalf)-(energiesShiftUpOneHalfLifetimeUnc);
%     
%     totalUncInShiftRelToScalar = sqrt(uncertaintyInShiftRelToScalarOneHalf.^2);
%     totalUncInMilliKelvinRelToScalar = prefactorTwo*totalUncInShiftRelToScalar/kB/1e-3;
%     totalUncInMHzRelToScalar = kB.*totalUncInMilliKelvinRelToScalar.*1e-3./hbar./2./pi./1e6;
%     
%     totalUncInShift = sqrt(uncertaintyInShiftTotalOneHalf.^2);
%     totalUncInMilliKelvin = prefactorTwo*totalUncInShift/kB/1e-3;
%     totalUncInMHz = kB.*totalUncInMilliKelvin.*1e-3./hbar./2./pi./1e6;
%     
%     totalEnergyShiftStateFactorShiftUpOneHalf = totalEnergyShiftStateFactor-normalEnergyOneHalf+energiesShiftUpOneHalfLifetimeUnc;
%     
%     alphaT = (totalEnergyShiftStateFactor(4)-totalEnergyShiftStateFactor(2))/(-totalEnergyShiftStateFactor(4));
%     alphaTUnc = (totalEnergyShiftStateFactorShiftUpOneHalf(4)-totalEnergyShiftStateFactorShiftUpOneHalf(2))/(-totalEnergyShiftStateFactorShiftUpOneHalf(4));
%     totalAlphaTUncertainty = sqrt((alphaTUnc-alphaT).^2);
%     
%     %     alphaVF1Down = (totalEnergyShiftStateFactor(3)-totalEnergyShiftStateFactor(1))/(2*totalEnergyShiftStateFactor(4));
%     alphaV = (totalEnergyShiftStateFactor(3)-totalEnergyShiftStateFactor(1))/(2*totalEnergyShiftStateFactor(4));
%     alphaVUnc = (totalEnergyShiftStateFactorShiftUpOneHalf(3)-totalEnergyShiftStateFactorShiftUpOneHalf(1))/(-2*totalEnergyShiftStateFactorShiftUpOneHalf(4));
%     totalAlphaVUncertainty = sqrt((alphaVUnc-alphaV).^2);
%     
%     figure(1232)
%     errorbar([-1:1:1],totalEnergyShiftInMHz(1:3)-totalEnergyShiftInMHz(4),totalUncInMHzRelToScalar(1:3),'.','MarkerSize',30);
%     xlim([-1.2 1.2])
%     
%     figure(1233)
%     errorbar([-1:1:1],-(totalEnergyShiftInMHz(1:3)-totalEnergyShiftInMHz(4))./totalEnergyShiftInMHz(4),totalUncInMHzRelToScalar(1:3)./totalEnergyShiftInMHz(4),'.','MarkerSize',30);
%     xlim([-1.2 1.2])
% end
