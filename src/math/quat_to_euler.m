function euler = quat_to_euler(q,p)
% Compute set od euler angles corresponding to the set of quaternions q
% according to the convention in p.sequence
% q can be 4 x N, euler = [ang3, ang2, ang1] shall be 3 x N
% The angles returned are ordered according to the convention
% RPY: [roll,pitch,yaw] (Default)
% ZYX: [yaw,pitch,roll]
% XYX (i.e. X1-Y-X2); [ang_x1,theta,ang_x2]
% REMARK: RPY and ZYX espress the same set of angles with the only difference
% of the order of the angles in the returned output

exp_str.sequence   = 'RPY';

if(~exist('p','var')), p = []; end
p = check_input_pars(exp_str,p);

if size(q,1)~=4
    error('Input must have 4 rows')
end

q0 = q(1,:);
q1 = q(2,:);
q2 = q(3,:);
q3 = q(4,:);

switch p.sequence
    case {'RPY','rpy'}
        
        sr_cp = 2*(q0.*q1 + q2.*q3);
        cr_cp = 1-2*(q1.^2 + q2.^2);
        sp = 2*(q0.*q2 - q3.*q1);
        sy_cp = 2*(q0.*q3 + q1.*q2);
        cy_cp = 1-2*(q2.^2 + q3.^2);
        
        r = atan2(sr_cp,cr_cp);        
        p = asin(sp);
        y = atan2(sy_cp,cy_cp);        
        
        euler = [r;p;y];
        
    case {'321','zyx','ZYX'}
        euler = dcm_to_euler(quat_to_dcm(q),p);
        
    case {'121','xyx','XYX'}
        euler = dcm_to_euler(quat_to_dcm(q),p);
        
    otherwise
        error(['Unknown sequence ' p.sequence]);
end

end