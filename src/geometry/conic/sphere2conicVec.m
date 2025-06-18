function conicVec = sphere2conicVec(pos_origin2sphere_REF, pos_sphere2cam_REF, dcm_REF2CAM, K, R)

conicMat = sphere2conicMat(pos_origin2sphere_REF, pos_sphere2cam_REF, dcm_REF2CAM, K, R);

conicVec = conicMat2conicVec(conicMat);

end