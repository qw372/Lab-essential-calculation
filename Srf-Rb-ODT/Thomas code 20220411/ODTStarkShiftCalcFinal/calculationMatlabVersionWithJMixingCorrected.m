%all energies in cm-1 unless otherwise specified.  All dipole moments in
%Debye unless otherwise specified

%important constants/conversions
hbar=1.0546e-34;
c=2.9979e8;
eps0=8.8542e-12;
debyeToMKS = 3.3356e-30;
invCmToHz=2.997e10;
e=1.602e-19;
kB=1.3806e-23;
bohrRadius = 5.2918e-11;
cmMinusOneToFrequency = 2.99793e10*2*pi;
debyeFactor = 0.3934;%convert debye to ea0

%'mystery factor'%
mysteryFactor = 1/2; %seemed to need a factor of 1/2 to make results match Varun+Yuqi. 
%At some point I think I reasoned out that this was necessary but now I'm not so sure...
%In any event, for now I'll leave this as a user adjustable parameter.  
%For what it's worth, the tensor/vector shift data (from microwave results)
%seem to either A) validate the 1/2 factor or B) imply that our intensity
%was actually 1/2 what we thought it was, given that they match nearly
%exactly with the 1/2 included, but removing it would make the tensor
%shifts also 2x greater...  

%flags
calculateUncertainties = 0;%calculate uncertainties in AC stark shift due to uncertainties in state lifetimes (and thus transition dipole moments), energy levels, etc.
addCoRotate=1; %we agreed that this is needed since 1064 nm is sufficiently far detuned such that negating co-rotating terms is not valid


%laser details
laserEnergy = 1/1064e-7;%for 1064 nm light
intensity = 2*51.5/(pi*(39.75e-6)^2);
p=1;%ODT polarization.  0 linear, +/-1 for \sigma^+\-

%molecule details

%J-mixing terms
a=0.888;
b=sqrt(1-a^2);

%electronic state energies
APiOneHalfEnergy = 15075.6122;
APiThreeHalvesEnergy = 15357.0736;
BSigmaEnergy = 17267.4465;
CPiEnergyOneHalf = 27384.67-57.9048/2;
CPiEnergyThreeHalves = 27384.67+57.9048/2;
DSigmaEnergy = 27773.8;
FSigmaEnergy = 32823.5;
GPiEnergy = 34808.9275;
HEnergy = 42345;

%detuning relative to each state energy, counter-rotating term
detuningAPiOneHalf = laserEnergy-APiOneHalfEnergy;
detuningAPiThreeHalves = laserEnergy-APiThreeHalvesEnergy;
detuningBSigma = laserEnergy-BSigmaEnergy;
detuningCPiOneHalf = laserEnergy-CPiEnergyOneHalf;
detuningCPiThreeHalves = laserEnergy-CPiEnergyThreeHalves;
detuningDSigma = laserEnergy-DSigmaEnergy;
detuningFSigma = laserEnergy-FSigmaEnergy;
detuningGPi = laserEnergy-GPiEnergy;
detuningH = laserEnergy-HEnergy;

%detuning relative to each state energy, co-rotating term
detuningAPiOneHalfCoRotate = laserEnergy+APiOneHalfEnergy;
detuningAPiThreeHalvesCoRotate = laserEnergy+APiThreeHalvesEnergy;
detuningBSigmaCoRotate = laserEnergy+BSigmaEnergy;
detuningCPiOneHalfCoRotate = laserEnergy+CPiEnergyOneHalf;
detuningCPiThreeHalvesCoRotate = laserEnergy+CPiEnergyThreeHalves;
detuningDSigmaCoRotate = laserEnergy+DSigmaEnergy;
detuningFSigmaCoRotate = laserEnergy+FSigmaEnergy;
detuningGPiCoRotate = laserEnergy+GPiEnergy;
detuningHCoRotate = laserEnergy+HEnergy;

%either calculate transition dipole moment based on experimentally measured
%lifetime data or from theoretical results for states for which there is no
%experimental data
debyeFromLifetime = @(lifetime,energy) sqrt(3*hbar/4/lifetime*4*pi*eps0*c^3/(energy*invCmToHz*2*pi)^3)/debyeToMKS;
debyeAPiOneHalf = debyeFromLifetime(24.1e-9+0e-9,APiOneHalfEnergy);%exp
debyeAPiThreeHalves = debyeFromLifetime(22.6e-9+0e-9,APiThreeHalvesEnergy);%exp
debyeBStates = debyeFromLifetime(25.5e-9-0.0e-9,BSigmaEnergy);%exp
debyeCPi = debyeFromLifetime(66.3e-9,CPiEnergyOneHalf);%de Melo+Ornellas
debyeDSigma = debyeFromLifetime(219e-9,DSigmaEnergy);%de Melo+Ornellas
debyeFSigma = 0.4/debyeFactor;%de Melo+Ornellas
debyeGPi = 0.5/debyeFactor;%de Melo + Ornellas
% debyeH = 2.2934;
debyeH = 0;

prefactor = mysteryFactor*(e^2*bohrRadius^2*debyeFactor^2)*intensity/2/c/eps0/hbar/cmMinusOneToFrequency;%First term is mystery factor.  Second term converts debye to MKS units.
%Remaining terms are necessary to convert (angular momenta matrix elements)^2/(detuning (in cm-1)) terms to trap depth.  

%X State values, in terms of |J,F,mF,N>
XValues = {[1/2,1,-1,1],[1/2,1,0,1],[1/2,1,1,1],[1/2,0,0,1],[3/2,1,-1,1],[3/2,1,0,1],[3/2,1,1,1],...
    [3/2,2,-2,1],[3/2,2,-1,1],[3/2,2,0,1],[3/2,2,1,1],[3/2,2,2,1]};

%%%FIRST do shifts due to APi_1/2 state.  This will involve calculations
%%%over J=1/2, J=3/2, and J=5/2 states |J,F,mF> in below

%J=1/2

APiOneHalfValues = {[1/2,1,-1],[1/2,1,0],[1/2,1,1],[1/2,0,0],...
    [3/2,1,-1],[3/2,1,0],[3/2,1,1],[3/2,2,-2],[3/2,2,-1],[3/2,2,0],[3/2,2,1],[3/2,2,2],...
    [5/2,2,-2],[5/2,2,-1],[5/2,2,0],[5/2,2,1],[5/2,2,2],[5/2,3,-3],...
    [5/2,3,-2],[5/2,3,-1],[5/2,3,0],[5/2,3,1],[5/2,3,2],[5/2,3,3]};

%APiOneHalfValues = {[1/2,1,-1],[1/2,1,0],[1/2,1,1],[1/2,0,0]};

%X state has
%(1/\sqrt(2))[-(1)^(j+1/2)|sigma=1/2,omega=1/2>+|sigma=-1/2,omega=-1/2>]
%APi_1/2 state has (1/\sqrt(2))[|sigma'=-1/2,omega'=1/2>+(-1)^(j'-1/2)|sigma'=-1/2,omega'=-1/2>]

XShiftDueToAPiOneHalf= zeros(1,12);

for i=1:length(XShiftDueToAPiOneHalf)
    Fp = XValues{i}(2);
    Jp = XValues{i}(1);
    for j=1:length(APiOneHalfValues)
        if Fp==0||Fp==2
            currAngularTermAbs = calculateAngularTerms(XValues{i},APiOneHalfValues{j},0,1/2,p);
            currAngularTermEmis = calculateAngularTerms(APiOneHalfValues{j},XValues{i},0,1/2,-p);
            XShiftDueToAPiOneHalf(i)=XShiftDueToAPiOneHalf(i)+abs(currAngularTermAbs)*abs(currAngularTermEmis);
            
            %special variation on calculation for dealing with angular
            %momentum mixing
        elseif Fp==1&&Jp==1/2
            currAngularTermMainState = calculateAngularTerms(XValues{i},APiOneHalfValues{j},0,1/2,p);
            currAngularTermMainStateEmis = calculateAngularTerms(APiOneHalfValues{j},XValues{i},0,1/2,-p);
            currAngularTermAdmixState = calculateAngularTerms([3/2,XValues{i}(2),XValues{i}(3)],APiOneHalfValues{j},0,1/2,p);
            XShiftDueToAPiOneHalf(i)=XShiftDueToAPiOneHalf(i)+(a*currAngularTermMainState-b*currAngularTermAdmixState)^2;
            
        elseif Fp==1&&Jp==3/2
            currAngularTermMainState = calculateAngularTerms(XValues{i},APiOneHalfValues{j},0,1/2,p);
            currAngularTermAdmixState = calculateAngularTerms([1/2,XValues{i}(2),XValues{i}(3)],APiOneHalfValues{j},0,1/2,p);
            XShiftDueToAPiOneHalf(i)=XShiftDueToAPiOneHalf(i)+(a*currAngularTermMainState+b*currAngularTermAdmixState)^2;
        end
    end
end

%%%NEXT do shifts due to A\Pi,Omega=3/2 state.  This will involve calculations
%%%over J=1/2, J=3/2, and J=5/2


APiThreeHalvesValues = {[1/2,1,-1],[1/2,1,0],[1/2,1,1],[1/2,0,0],...
    [3/2,1,-1],[3/2,1,0],[3/2,1,1],[3/2,2,-2],[3/2,2,-1],[3/2,2,0],[3/2,2,1],[3/2,2,2],...
    [5/2,2,-2],[5/2,2,-1],[5/2,2,0],[5/2,2,1],[5/2,2,2],[5/2,3,-3],...
    [5/2,3,-2],[5/2,3,-1],[5/2,3,0],[5/2,3,1],[5/2,3,2],[5/2,3,3]};

%X state has
%(1/\sqrt(2))[-(1)^(j+1/2)|sigma=1/2,omega=1/2>+|sigma=-1/2,omega=-1/2>]
%APi_1/2 state has (1/\sqrt(2))[|sigma'=-1/2,omega'=1/2>+(-1)^(j'-1/2)|sigma'=-1/2,omega'=-1/2>]

XShiftDueToAPiThreeHalves = zeros(1,12);

for i=1:length(XShiftDueToAPiThreeHalves)
    Fp = XValues{i}(2);
    Jp = XValues{i}(1);
    for j=1:length(APiThreeHalvesValues)
        if Fp==0||Fp==2
            currAngularTermAbs = calculateAngularTerms(XValues{i},APiThreeHalvesValues{j},0,3/2,p);
            XShiftDueToAPiThreeHalves(i)=XShiftDueToAPiThreeHalves(i)+currAngularTermAbs^2;
            
        elseif Fp==1&&Jp==1/2
            currAngularTermMainState = calculateAngularTerms(XValues{i},APiThreeHalvesValues{j},0,3/2,p);
            currAngularTermAdmixState = calculateAngularTerms([3/2,XValues{i}(2),XValues{i}(3)],APiThreeHalvesValues{j},0,3/2,p);
            XShiftDueToAPiThreeHalves(i)=XShiftDueToAPiThreeHalves(i)+(a*currAngularTermMainState-b*currAngularTermAdmixState)^2;
            
        elseif Fp==1&&Jp==3/2
            currAngularTermMainState = calculateAngularTerms(XValues{i},APiThreeHalvesValues{j},0,3/2,p);
            currAngularTermAdmixState = calculateAngularTerms([1/2,XValues{i}(2),XValues{i}(3)],APiThreeHalvesValues{j},0,3/2,p);
            XShiftDueToAPiThreeHalves(i)=XShiftDueToAPiThreeHalves(i)+(a*currAngularTermMainState+b*currAngularTermAdmixState)^2;
        end
    end
end

%%%NEXT do shifts due to B\Sigma state.  This will involve calculations
%%%over N=0 and N=2 states

%N=0 (B sigma ordering: [J,F,mF,N]

BSigmaValues = {[1/2,1,-1,0],[1/2,1,0,0],[1/2,1,1,0],[1/2,0,0,0],...
    [3/2,1,-1,2],[3/2,1,0,2],[3/2,1,1,2],[3/2,2,-2,2],[3/2,2,-1,2],[3/2,2,0,2],...
    [3/2,2,1,2],[3/2,2,2,2],[5/2,2,-2,2],[5/2,2,-1,2],[5/2,2,0,2],[5/2,2,1,2],[5/2,2,2,2],...
    [5/2,3,-3,2],[5/2,3,-2,2],[5/2,3,-1,2],[5/2,3,0,2],[5/2,3,1,2],[5/2,3,2,2],[5/2,3,3,2]};

XShiftDueToBSigma = zeros(1,12);

for i=1:length(XShiftDueToBSigma)
    Fp = XValues{i}(2);
    Jp = XValues{i}(1);
    for j=1:length(BSigmaValues)
        if Fp==0||Fp==2
            currAngularTermAbs = calculateAngularTerms(XValues{i},BSigmaValues{j},1,1/2,p);
            XShiftDueToBSigma(i)=XShiftDueToBSigma(i)+currAngularTermAbs^2;
            
        elseif Fp==1&&Jp==1/2
            currAngularTermMainState = calculateAngularTerms(XValues{i},BSigmaValues{j},1,1/2,p);
            currAngularTermAdmixState = calculateAngularTerms([3/2,XValues{i}(2),XValues{i}(3),XValues{i}(4)],BSigmaValues{j},1,1/2,p);
            XShiftDueToBSigma(i)=XShiftDueToBSigma(i)+(a*currAngularTermMainState-b*currAngularTermAdmixState)^2;
            
        elseif Fp==1&&Jp==3/2
            currAngularTermMainState = calculateAngularTerms(XValues{i},BSigmaValues{j},1,1/2,p);
            currAngularTermAdmixState = calculateAngularTerms([1/2,XValues{i}(2),XValues{i}(3),XValues{i}(4)],BSigmaValues{j},1,1/2,p);
            XShiftDueToBSigma(i)=XShiftDueToBSigma(i)+(a*currAngularTermMainState+b*currAngularTermAdmixState)^2;
        end
    end
end

totalEnergyShiftStateFactor = (debyeAPiOneHalf^2/detuningAPiOneHalf).*(XShiftDueToAPiOneHalf)+...
    (debyeAPiThreeHalves^2/detuningAPiThreeHalves).*(XShiftDueToAPiThreeHalves)+...
    (debyeBStates^2/detuningBSigma).*(XShiftDueToBSigma)+...
    (debyeCPi^2/detuningCPiOneHalf).*(XShiftDueToAPiOneHalf)+...%angular momenta calcs are same for all Pi states and all Sigma states, so only need to do it once for 'A' and 'B' respectively and then apply to other states
    (debyeCPi^2/detuningCPiThreeHalves).*(XShiftDueToAPiThreeHalves)+...
    (debyeDSigma^2/detuningDSigma).*(XShiftDueToBSigma)+...
    (debyeFSigma^2/detuningFSigma).*(XShiftDueToBSigma)+...
    (debyeGPi^2/detuningGPi).*(XShiftDueToAPiOneHalf+XShiftDueToAPiThreeHalves)+...
    (debyeH^2/detuningH)*1/3;

if addCoRotate==1%NOTE: minus signs are correct here, see ODT grimm paper...really the counter-rot detuning should be pos and the signs above should also be negative, see Grimm ODT review
    totalEnergyShiftStateFactor = -(debyeAPiOneHalf^2.*(1./(-laserEnergy+APiOneHalfEnergy)+1./(laserEnergy+APiOneHalfEnergy))).*(XShiftDueToAPiOneHalf)-...
        (debyeAPiThreeHalves^2.*(1./(-laserEnergy+APiThreeHalvesEnergy)+1./(laserEnergy+APiThreeHalvesEnergy))).*(XShiftDueToAPiThreeHalves)-...
        (debyeBStates^2.*(1./(-laserEnergy+BSigmaEnergy)+1./(laserEnergy+BSigmaEnergy))).*(XShiftDueToBSigma)-...
        (debyeCPi^2.*(1./(-laserEnergy+CPiEnergyOneHalf)+1./(laserEnergy+CPiEnergyOneHalf))).*(XShiftDueToAPiOneHalf)-...
        (debyeCPi^2.*(1./(-laserEnergy+CPiEnergyThreeHalves)+1./(laserEnergy+CPiEnergyThreeHalves))).*(XShiftDueToAPiThreeHalves)-...
        (debyeDSigma^2.*(1./(-laserEnergy+DSigmaEnergy)+1./(laserEnergy+DSigmaEnergy))).*(XShiftDueToBSigma)-...
        (debyeFSigma^2.*(1./(-laserEnergy+FSigmaEnergy)+1./(laserEnergy+FSigmaEnergy))).*(XShiftDueToBSigma)-...
        (debyeGPi^2.*(1./(-laserEnergy+GPiEnergy)+1./(laserEnergy+GPiEnergy))).*(XShiftDueToAPiOneHalf+XShiftDueToAPiThreeHalves)-...
        (debyeH^2.*(1./(-laserEnergy+HEnergy)+1./(laserEnergy+HEnergy)))*1/3;
end


%results
%Trap depth
totalEnergyShiftInMilliKelvin2 = (prefactor)*totalEnergyShiftStateFactor/kB/1e-3
totalEnergyShiftInMHz = kB.*totalEnergyShiftInMilliKelvin2.*1e-3./hbar./2./pi./1e6;%total energy shifts of all 12 |X\Sigma> hyperfine states
%Tensor/Vector shifts
alphaTF1Down = (totalEnergyShiftStateFactor(4)-totalEnergyShiftStateFactor(2))/(-totalEnergyShiftStateFactor(4));
alphaVF1Down = (totalEnergyShiftStateFactor(3)-totalEnergyShiftStateFactor(1))/(2*totalEnergyShiftStateFactor(4));
alphaTF1Up = (totalEnergyShiftStateFactor(4)-totalEnergyShiftStateFactor(6))/(-totalEnergyShiftStateFactor(4));
alphaVF1Up = (totalEnergyShiftStateFactor(7)-totalEnergyShiftStateFactor(5))/(2*totalEnergyShiftStateFactor(4));
alphaTF2 = 2*(totalEnergyShiftStateFactor(4)-totalEnergyShiftStateFactor(10))/(-totalEnergyShiftStateFactor(4));
alphaVF2 = (totalEnergyShiftStateFactor(12)-totalEnergyShiftStateFactor(8))/(2*totalEnergyShiftStateFactor(4));
    
if calculateUncertainties ==1
    aPiOneHalfLifetimeUnc = -2/24.1*1;
    aPiThreeHalfLifetimeUnc = -4.7/22.6*1;
    
    energiesShiftUpOneHalfLifetimeUnc =  (debyeAPiOneHalf^2/detuningAPiOneHalf).*(XShiftDueToAPiOneHalf)*(1+aPiOneHalfLifetimeUnc);
    normalEnergyOneHalf = (debyeAPiOneHalf^2/detuningAPiOneHalf).*(XShiftDueToAPiOneHalf)*(1);
    uncertaintyInShiftRelToScalarOneHalf = (normalEnergyOneHalf-normalEnergyOneHalf(4))-(energiesShiftUpOneHalfLifetimeUnc-energiesShiftUpOneHalfLifetimeUnc(4));
    uncertaintyInShiftTotalOneHalf = (normalEnergyOneHalf)-(energiesShiftUpOneHalfLifetimeUnc);
    
    energiesShiftUpThreeHalfLifetimeUnc =  (debyeAPiThreeHalves^2/detuningAPiThreeHalves).*(XShiftDueToAPiThreeHalves)*(1+aPiThreeHalfLifetimeUnc);
    normalEnergyThreeHalf = (debyeAPiThreeHalves^2/detuningAPiThreeHalves).*(XShiftDueToAPiThreeHalves)*(1);
    uncertaintyInShiftRelToScalarThreeHalf = (normalEnergyThreeHalf-normalEnergyThreeHalf(4))-(energiesShiftUpThreeHalfLifetimeUnc-energiesShiftUpThreeHalfLifetimeUnc(4));
    uncertaintyInShiftTotalThreeHalf = (normalEnergyThreeHalf)-(energiesShiftUpThreeHalfLifetimeUnc);
    
    totalUncInShiftRelToScalar = sqrt(uncertaintyInShiftRelToScalarOneHalf.^2+uncertaintyInShiftRelToScalarThreeHalf.^2);
    totalUncInMilliKelvinRelToScalar = (prefactor)*totalUncInShiftRelToScalar/kB/1e-3;
    totalUncInMHzRelToScalar = kB.*totalUncInMilliKelvinRelToScalar.*1e-3./hbar./2./pi./1e6;
    
    totalUncInShift = sqrt(uncertaintyInShiftTotalOneHalf.^2+uncertaintyInShiftTotalThreeHalf.^2);
    totalUncInMilliKelvin = (prefactor)*totalUncInShift/kB/1e-3;
    totalUncInMHz = kB.*totalUncInMilliKelvin.*1e-3./hbar./2./pi./1e6;
    
    totalEnergyShiftStateFactorShiftUpOneHalf = totalEnergyShiftStateFactor-normalEnergyOneHalf+energiesShiftUpOneHalfLifetimeUnc;
    totalEnergyShiftStateFactorShiftUpThreeHalf = totalEnergyShiftStateFactor-normalEnergyThreeHalf+energiesShiftUpThreeHalfLifetimeUnc;
    
    
    alphaTF1DownOneHalfUpUnc = (totalEnergyShiftStateFactorShiftUpOneHalf(4)-totalEnergyShiftStateFactorShiftUpOneHalf(2))/(-totalEnergyShiftStateFactorShiftUpOneHalf(4));
    alphaTF1DownThreeHalfUpUnc = (totalEnergyShiftStateFactorShiftUpThreeHalf(4)-totalEnergyShiftStateFactorShiftUpThreeHalf(2))/(-totalEnergyShiftStateFactorShiftUpThreeHalf(4));
    totalAlphaTF1DownUncertainty = sqrt((alphaTF1DownOneHalfUpUnc-alphaTF1Down).^2+(alphaTF1DownThreeHalfUpUnc-alphaTF1Down).^2);
    
    
    alphaVF1DownOneHalfUpUnc = (totalEnergyShiftStateFactorShiftUpOneHalf(3)-totalEnergyShiftStateFactorShiftUpOneHalf(1))/(-2*totalEnergyShiftStateFactorShiftUpOneHalf(4));
    alphaVF1DownThreeHalfUpUnc = (totalEnergyShiftStateFactorShiftUpThreeHalf(3)-totalEnergyShiftStateFactorShiftUpThreeHalf(1))/(-2*totalEnergyShiftStateFactorShiftUpThreeHalf(4));
    totalAlphaVF1DownUncertainty = sqrt((alphaVF1DownOneHalfUpUnc-alphaVF1Down).^2+(alphaVF1DownThreeHalfUpUnc-alphaVF1Down).^2);
    
    
    alphaTF1UpOneHalfUpUnc = (totalEnergyShiftStateFactorShiftUpOneHalf(4)-totalEnergyShiftStateFactorShiftUpOneHalf(6))/(-totalEnergyShiftStateFactorShiftUpOneHalf(4));
    alphaTF1UpThreeHalfUpUnc = (totalEnergyShiftStateFactorShiftUpThreeHalf(4)-totalEnergyShiftStateFactorShiftUpThreeHalf(6))/(-totalEnergyShiftStateFactorShiftUpThreeHalf(4));
    totalAlphaTF1UpUncertainty = sqrt((alphaTF1UpOneHalfUpUnc-alphaTF1Up).^2+(alphaTF1UpThreeHalfUpUnc-alphaTF1Up).^2);
    
    
    alphaVF1UpOneHalfUpUnc = (totalEnergyShiftStateFactorShiftUpOneHalf(7)-totalEnergyShiftStateFactorShiftUpOneHalf(5))/(-2*totalEnergyShiftStateFactorShiftUpOneHalf(4));
    alphaVF1UpThreeHalfUpUnc = (totalEnergyShiftStateFactorShiftUpThreeHalf(7)-totalEnergyShiftStateFactorShiftUpThreeHalf(5))/(-2*totalEnergyShiftStateFactorShiftUpThreeHalf(4));
    totalAlphaVF1UpUncertainty = sqrt((alphaVF1UpOneHalfUpUnc-alphaVF1Up).^2+(alphaVF1UpThreeHalfUpUnc-alphaVF1Up).^2);
    
    
    
    figure(1232)
    errorbar([-1:1:1],totalEnergyShiftInMHz(1:3)-totalEnergyShiftInMHz(4),totalUncInMHzRelToScalar(1:3),'.','MarkerSize',30);
    hold all
    errorbar([-1:1:1],totalEnergyShiftInMHz(5:7)-totalEnergyShiftInMHz(4),totalUncInMHzRelToScalar(5:7),'.','MarkerSize',30);
    errorbar([-2:1:2],totalEnergyShiftInMHz(8:12)-totalEnergyShiftInMHz(4),totalUncInMHzRelToScalar(8:12),'.','MarkerSize',30);
    xlim([-2.2 2.2])
    
    figure(1233)
    errorbar([-1:1:1],-(totalEnergyShiftInMHz(1:3)-totalEnergyShiftInMHz(4))./totalEnergyShiftInMHz(4),totalUncInMHzRelToScalar(1:3)./totalEnergyShiftInMHz(4),'.','MarkerSize',30);
    hold all
    errorbar([-1:1:1],-(totalEnergyShiftInMHz(5:7)-totalEnergyShiftInMHz(4))./totalEnergyShiftInMHz(4),totalUncInMHzRelToScalar(5:7)./totalEnergyShiftInMHz(4),'.','MarkerSize',30);
    errorbar([-2:1:2],-(totalEnergyShiftInMHz(8:12)-totalEnergyShiftInMHz(4))./totalEnergyShiftInMHz(4),totalUncInMHzRelToScalar(8:12)./totalEnergyShiftInMHz(4),'.','MarkerSize',30);
    xlim([-2.2 2.2])
end
