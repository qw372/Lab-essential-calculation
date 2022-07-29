function angularTerm = calculateAngularTermsAPi(GndVals,ExcVals,ExcAOrB,ExcOmega,p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%assumes ground state omega is 1/2 and N=1, would need to rewrite last bit if you
%want to make it more adaptable to arbitrary omega

if ExcAOrB==0 %0 if going from A to a case (b) state
    
    currAVals = GndVals;
    Jp=currAVals(1);
    Fp=currAVals(2);
    mFp=currAVals(3);
    currDelVals = ExcVals;
    J = currDelVals(1);
    F = currDelVals(2);
    mF = currDelVals(3);
    termOne = sqrt((2*F+1)*(2*Fp+1)*(2*J+1)*(2*Jp+1));
    termTwo = (-1)^(F-mF)*Wigner3j([F,1,Fp],[-mF,p,mFp]);
    if (F==0&&Fp==0 || abs(F-Fp)>1 || abs(J-Jp)>1)
        termThree = 0;
    else
        termThree = (-1)^(Fp+J+1+1/2)*Wigner6j(Jp,Fp,1/2,F,J,1);
    end
    if ExcOmega==1/2
        termFour = (1/2*((-1)^(J+Jp)*(-1)^(J+1/2)*Wigner3j([J,1,Jp],[1/2,-1,1/2])+(-1)^(J-1/2)*Wigner3j([J,1,Jp],[-1/2,1,-1/2])));
    elseif ExcOmega==3/2
        termFour = (1/2*((-1)^(Jp+1/2)*(-1)^(J-3/2)*Wigner3j([J,1,Jp],[-3/2,1,1/2])+(-1)^(J-1/2)*(-1)^(J+3/2)*Wigner3j([J,1,Jp],[3/2,-1,-1/2])));
    end
    angularTerm = termOne*termTwo*termThree*termFour;
    
elseif ExcAOrB==1 %this is an Pi to Pi transition in this case
    
    currAVals = GndVals;
    Jp=currAVals(1);
    Fp=currAVals(2);
    mFp=currAVals(3);
    currCVals = ExcVals;
    J = currCVals(1);
    F = currCVals(2);
    mF = currCVals(3);
    termOne = sqrt((2*F+1)*(2*Fp+1)*(2*J+1)*(2*Jp+1));
    termTwo = (-1)^(F-mF)*Wigner3j([F,1,Fp],[-mF,p,mFp]);
    if (F==0&&Fp==0 || abs(F-Fp)>1 || abs(J-Jp)>1)
        termThree = 0;
    else
        termThree = (-1)^(Fp+J+1+1/2)*Wigner6j(Jp,Fp,1/2,F,J,1);
    end
    termFour = (1/2*((-1)^(J-1/2)*Wigner3j([J,1,Jp],[-1/2,0,1/2])+(-1)^(J+1/2)*(-1)^(J+1-1/2)*Wigner3j([J,1,Jp],[1/2,0,-1/2])));
    
    
    angularTerm = termOne*termTwo*termThree*termFour;
    
    
elseif ExcAOrB==2 %this is an Pi to Delta transition in this case
    
    currAVals = GndVals;
    Jp=currAVals(1);
    Fp=currAVals(2);
    mFp=currAVals(3);
    currDelVals = ExcVals;
    J = currDelVals(1);
    F = currDelVals(2);
    mF = currDelVals(3);
    termOne = sqrt((2*F+1)*(2*Fp+1)*(2*J+1)*(2*Jp+1));
    termTwo = (-1)^(F-mF)*Wigner3j([F,1,Fp],[-mF,p,mFp]);
    if (F==0&&Fp==0 || abs(F-Fp)>1 || abs(J-Jp)>1)
        termThree = 0;
    else
        termThree = (-1)^(Fp+J+1+1/2)*Wigner6j(Jp,Fp,1/2,F,J,1);
    end
    termFour = (1/2*((-1)^(J-3/2)*Wigner3j([J,1,Jp],[-3/2,1,1/2])+(-1)^(J+3/2)*(-1)^(J+1-1/2)*Wigner3j([J,1,Jp],[3/2,-1,-1/2])));
    
    
    angularTerm = termOne*termTwo*termThree*termFour;
    

end

