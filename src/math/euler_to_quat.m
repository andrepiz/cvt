function quat = euler_to_quat(euler, p)
% EULER_TO_QUAT This function converts from euler angles to a quaternions 
% Input: Euler angles [3 x n-set-of-angles] according to p.sequence convention [rad]
% Outputs: quaternions [4 x n-set-of-angles]
% The angles in input must be ordered according to the convention
% RPY: [roll,pitch,yaw] (Default)
% ZYX: [yaw,pitch,roll]
% REMARK: RPY and ZYX espress the same set of angles with the only difference
% of the assumed order of the angles in the input

exp_str.sequence   = 'RPY';

if(~exist('p','var')), p = []; end
p = check_input_pars(exp_str,p);

if size(euler,1)~=3
    error('Input must have 3 rows')
end

switch p.sequence
    case {'RPY','rpy'}
        s = sin(0.5*euler);
        c = cos(0.5*euler);
        z0 = zeros(size(euler(1,:,:)));
        q1 = [c(3,:,:); z0; z0; s(3,:,:)];
        q2 = [c(2,:,:); z0; s(2,:,:); z0];
        q3 = [c(1,:,:); s(1,:,:); z0; z0];
        quat = quat_mult(quat_mult(q1,q2),q3);
    case {'ZYX','zyx'}
        quat = euler_to_quat(euler([3,2,1]',:,:));
    case {'YZX','yzx'}
        quat = euler_to_quat(euler([2,3,1]',:,:));
    otherwise
        quat = dcm_to_quat(euler_to_dcm(euler,p));
end

end
    