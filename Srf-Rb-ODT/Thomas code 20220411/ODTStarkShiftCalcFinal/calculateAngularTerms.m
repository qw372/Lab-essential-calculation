function angularTerm = calculateAngularTerms(GndVals,ExcVals,ExcAOrB,ExcOmega,p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%assumes ground state omega is 1/2 and N=1, would need to rewrite last bit if you
%want to make it more adaptable to arbitrary omega

if ExcAOrB==0
    
    currXVals = GndVals;
    Jp=currXVals(1);
    Fp=currXVals(2);
    mFp=currXVals(3);
    currAVals = ExcVals;
    J = currAVals(1);
    F = currAVals(2);
    mF = currAVals(3);
    if length(currXVals)==4
        flip=0;
        JX=Jp;
        JA=J;
    else
        flip=1;
        JX=J;
        JA=Jp;
    end
    termOne = sqrt((2*F+1)*(2*Fp+1)*(2*J+1)*(2*Jp+1));
    termTwo = (-1)^(F-mF)*Wigner3j([F,1,Fp],[-mF,p,mFp]);
    if (F==0&&Fp==0 || abs(F-Fp)>1 || abs(J-Jp)>1)
        termThree = 0;
    else
        termThree = (-1)^(Fp+J+1+1/2)*Wigner6j(Jp,Fp,1/2,F,J,1);
    end
    if ExcOmega==1/2
          termFour = (1/2*((-1)^(J+Jp)*(-1)^(J+1/2)*Wigner3j([J,1,Jp],[1/2,-1,1/2])+(-1)^(J-1/2)*Wigner3j([J,1,Jp],[-1/2,1,-1/2])));
%         if flip==0
%             termFour = 1/2*((-1)^(JX+1/2+JA-1/2+J-1/2)*Wigner3j([J,1,Jp],[1/2,-1,1/2])+(-1)^(J+1/2)*Wigner3j([J,1,Jp],[-1/2,1,-1/2]));
%         else
%             termFour = 1/2*((-1)^(J-1/2)*Wigner3j([J,1,Jp],[1/2,-1,1/2])+(-1)^(JX+1/2+JA-1/2+J+1/2)*Wigner3j([J,1,Jp],[-1/2,1,-1/2]));
%         end
    elseif ExcOmega==3/2
        termFour = (1/2*((-1)^(Jp+1/2)*(-1)^(J-3/2)*Wigner3j([J,1,Jp],[-3/2,1,1/2])+(-1)^(J-1/2)*(-1)^(J+3/2)*Wigner3j([J,1,Jp],[3/2,-1,-1/2])));
    end
    angularTerm = termOne*termTwo*termThree*termFour;
    
elseif ExcAOrB==1
    
    currXVals = GndVals;
    Jp=currXVals(1);
    Fp=currXVals(2);
    mFp=currXVals(3);
    Np = currXVals(4);
    currBVals = ExcVals;
    J = currBVals(1);
    F = currBVals(2);
    mF = currBVals(3);
    N = currBVals(4);
    
    termOne = sqrt((2*F+1)*(2*Fp+1)*(2*J+1)*(2*Jp+1)*(2*Np+1)*(2*N+1));
    termTwo = (-1)^(F-mF)*Wigner3j([F,1,Fp],[-mF,p,mFp]);
    
    if (F==0&&Fp==0 || abs(F-Fp)>1 || abs(J-Jp)>1)
        termThree = 0;
    else
        termThree = (-1)^(Fp+J+1+1/2)*Wigner6j(Jp,Fp,1/2,F,J,1);
    end
    
    if (abs(J-Jp)>1)
        termFour = 0;
    else
        termFour = (-1)^(Jp+N+1+1/2)*Wigner6j(Np,Jp,1/2,J,N,1);
    end
    
    termFive = (-1)^N*Wigner3j([N,1,Np],[0,0,0]);
    angularTerm = termOne*termTwo*termThree*termFour*termFive;
end

