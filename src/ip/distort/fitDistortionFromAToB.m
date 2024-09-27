function [coeff,tab]=fitDistortionFromAToB(LOS_A,LOS_B,polyDegree,effIrr)
LOS_A = reshape(LOS_A,3,[],1); 
LOS_B = reshape(LOS_B,3,[],1); 

LOS_A_h= LOS_A(1:2,:)./LOS_A(3,:);
LOS_B_h= LOS_B(1:2,:)./LOS_B(3,:);

idxOk= all(~isnan(LOS_A_h)) & all(~isnan(LOS_B_h));
if nargin<4
useIrrCorrection=false;
effIrr =1; 
else 
useIrrCorrection=true;
effIrr = reshape(effIrr,1,[],1);
effIrr = effIrr(idxOk);
end

LOS_A_h = LOS_A_h(:,idxOk);
LOS_B_h=LOS_B_h(:,idxOk);
nSamples = size(LOS_A_h,2);

% build degrees table
[nCoeff,tab]= assembleTab(polyDegree,useIrrCorrection);

% build least squares matrix A
A = zeros(nSamples,nCoeff);
for kk=1:nCoeff
    A(:,kk)=(LOS_A_h(1,:).^tab(kk,1)).*(LOS_A_h(2,:).^tab(kk,2)).*(effIrr.^tab(kk,3)); 
end

% get solution
%solS=;
coeff=((A'*A)\A')*LOS_B_h';
 
end



function [numCoeff,tab]= assembleTab(polyDegree,useIrradianceCorrection)


numCoeff=0;
for ii=polyDegree:-1:0
    for jj=polyDegree-ii:-1:0
        if useIrradianceCorrection
            for kk=polyDegree-ii-jj:-1:0
                numCoeff=numCoeff+1;
            end
        else
            numCoeff=numCoeff+1;
        end
    end
end


tab=zeros(numCoeff,3);
mm=0;
for ii=polyDegree:-1:0
    for jj=polyDegree-ii:-1:0
        if useIrradianceCorrection
            for kk=polyDegree-ii-jj:-1:0
                mm=mm+1;
                tab(mm,:)=[ii,jj,kk];
            end
        else
            mm=mm+1;
            tab(mm,:)=[ii,jj,0];
        end
    end
end


end