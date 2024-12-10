function euler = dcm_to_euler(M,p)
% Compute set od euler angles corresponding to the set of DCM matrices M
% according to the convention in p.sequence
% M can be 3 x 3 or 3 x 3 x N, euler = [ang1, ang2, ang3] shall be 3 x N
% The angles returned are ordered according to the convention
% RPY: [roll,pitch,yaw] (Default)
% ZYX: [yaw,pitch,roll]
% XYX (i.e. X1-Y-X2); [ang_x1,theta,ang_x2]
% REMARK: RPY and ZYX espress the same set of angles with the only difference
% of the order of the angles in the returned output

exp_str.sequence   = 'RPY';

if(~exist('p','var')), p = []; end
p = check_input_pars(exp_str,p);

if size(M,1)~=3 || size(M,2)~=3
    error('Input must be 3 x 3 x N')
end

switch p.sequence
    case {'RPY','rpy'}
        m11 = M(1,1,:);
        m12 = M(1,2,:);
        m13 = M(1,3,:);
        m23 = M(2,3,:);
        m33 = M(3,3,:);
        % Finding phi
        roll = atan2(m23,m33);
        % Finding psi
        yaw = atan2(m12,m11);
        % Finding theta
        pitch = asin(-m13);
        
        euler = [roll(:),pitch(:),yaw(:)]';        
    case {'321','zyx','ZYX'}
        m11 = M(1,1,:);
        m12 = M(1,2,:);
        m13 = M(1,3,:);
        m23 = M(2,3,:);
        m33 = M(3,3,:);
        % Finding phi
        roll = atan2(m23,m33);
        % Finding psi
        yaw = atan2(m12,m11);
        % Finding theta
        pitch = asin(-m13);
        
        euler = [yaw(:),pitch(:),roll(:)]';
    case {'121','xyx','XYX'}
        m11 = M(1,1,:);
        m12 = M(1,2,:);
        m13 = M(1,3,:);
        m21 = M(2,1,:);
        m31 = M(3,1,:);
        % Finding phi
        phi = atan2(m12,-m13);
        % Finding psi
        psi = atan2(m21,m31);
        % Finding theta
        theta = acos(m11);
        
        euler = [phi(:),theta(:),psi(:)]';
    otherwise
        error(['Unknown sequence ' p.sequence]);
end

end