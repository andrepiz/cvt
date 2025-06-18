function conicVec = conicMat2conicVec(conicMat)

conicVec(1) = conicMat(1,1);
conicVec(2) = 2*conicMat(1,2);
conicVec(3) = conicMat(2,2);
conicVec(4) = 2*conicMat(1,3);
conicVec(5) = 2*conicMat(2,3);
conicVec(6) = conicMat(3,3);

end

