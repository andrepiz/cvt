function x_rot_in_A = rotvector(x_in_A, q_rot, conj)
% ROTVECTOR Assumes x vector components are in A frame and returns the components in
% that same frame if x is rotated according to the quaternion q 
% This is equivalent to the convention fixed-frame-rotating-vector.
% If conj is true, the quaternion is conjugated so the rotations occurs
% in the contrary direction as specified by q.

if(~exist('conj','var'))
    conj = false;
end

if ~(size(q_rot,1)==4)
    error('quaternion is not 4-by-N')
end
if ~(size(x_in_A,1)==3)
    error('vector is not 3-by-N')
end
N1 = size(x_in_A,2);
N2 = size(q_rot,2);

if N1~=N2 && N1~=1 && N2~=1
    error('Vector and quaternion sizes are not compatible.')
end

if(conj)
    q_rot = quat_conj(q_rot);
end

if N1==1
    x_in_A = x_in_A*ones(1,N2);
end
if N2==1
    q_rot = q_rot*ones(1,N1);
end

q0 = q_rot(1,:);
qv = q_rot(2:end,:);
x_rot_in_A = x_in_A + cross(2*qv, cross(qv, x_in_A) + repmat(q0,3,1).*x_in_A);

end