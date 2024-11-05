function x_B = rotframe(x_A, q_AtoB, conj)
% ROTFRAME Assumes x vector components are in A frame and returns the
% components in B frame using the quaternion rotation q from A
% to B. This is equivalent to the convention fixed-vector-rotating-frame.
% If conj is true, the quaternion is conjugated so the method
% assumes the components are in B frame and are converted to A frame.
%
% x_B = rotframe(x_A,q_AtoB,conj) express x_A vector [3 x N1] expressed in frame A
%                                    in the equivalent x_B vector in frame B whre  q_AtoB [4 x N2]
%                                    is the rotation from A to B
% possibilities N1=N2, N1=1 & N2 whatever,  N1 qhatever & N2=1

if ~(size(q_AtoB,1)==4)
    error('quaternion is not 4-by-N')
end
if ~(size(x_A,1)==3)
    error('vector is not 3-by-N')
end

N1 = size(x_A,2);
N2 = size(q_AtoB,2);

if N1~=N2 && N1~=1 && N2~=1
    error('Vector and quaternion sizes are not compatible.')
end


if(~exist('conj','var'))
    conj = false;
end

if(conj)
    % If conj is true, x_A is actually interpreted as x_B. q_AtoB is the quaternion from A to B;
    % we want the rotation to be from B to A, 
    % so it is needed to calculate the q_BtoA (as the conj of q_AtoB). 
    q_BtoA = quat_conj(q_AtoB);
    x_B = x_A;
    % the operation is equivalent to the fixed-frame-rotating-vector operation
    % but with conjugated quaternion
    x_A = rotframe(x_B, q_BtoA, false);
    x_B = x_A;
else
    x_B = rotvector(x_A, q_AtoB, true);
end



end