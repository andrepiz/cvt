function [ix_blind_fov, ix_blind_baffle, sc2body_az_urf, sc2body_ze_urf, sepAngReq_fov, sepAngReq_baffle] = ...
    check_visibility(this, bodies_list, et, rr_eci, q_eci2rbf, p)
    % Checks if point observed from a camera is in the FoV, comparing
    % the angle of the celestial body with respect to the boresight 
    % with the sum of the angular dimension of the sphere
    % and the field of view of the sensor.
    % Output are in [nb x np x nd] dimension, where:
    %   nb: number of celestial bodies
    %   np: number of points
    %   nd: number of STRs

    exp_str.fov_shape           = 'circle';
    exp_str.baffle_shape        = 'circle';
    exp_str.observer            = 'earth';
    exp_str.rpy_offpoint        = [0;0;0];

    if(~exist('p','var')), p = []; end
    p = cp_check_input_pars(exp_str,p);
    
    nb = length(bodies_list);
    nd = length(this);
    np = size(rr_eci, 2);

    % init
    sepAngReq_fov = zeros(nb, np, nd);
    sepAngReq_baffle = zeros(nb, np, nd);
    sc2body_az_urf = zeros(nb, np, nd);
    sc2body_ze_urf = zeros(nb, np, nd);
    ix_blind_fov = zeros(nb, np, nd);
    ix_blind_baffle = zeros(nb, np, nd);
    obs2body_eci = zeros(3,np,nb);
    sc2body_eci  = zeros(3,np,nb);
    sc2body_norm_eci = zeros(1,np,nb);
    sc2body_dir_eci = zeros(3,np,nb);
    R_body = zeros(1,np);
    angDim_body  = zeros(1,np,nb);
    q_rbf2urf = zeros(4,nd);
    q_eci2urf = zeros(4,np,nd);
    sc2body_dir_urf  = zeros(3,np,nb,nd);
    sc2body_sph_urf  = zeros(3,np,nb,nd);
    
    for i = 1:nb

        % direction and semicone angle of each body wrt sc
        obs2body_eci(:,:,i) = aca_body_position(bodies_list{i}, p.observer, 'J2000', et);
        sc2body_eci(:,:,i) = -rr_eci + obs2body_eci(:,:,i);
        sc2body_norm_eci(1,:,i) = cp_norm(sc2body_eci(:,:,i));
        sc2body_dir_eci(:,:,i) = cp_normalize(sc2body_eci(:,:,i));
        R_body(i) = Environment.Body.(bodies_list{i}).R_v;
        angDim_body(1,:,i) = asin(R_body(i)./sc2body_norm_eci(1,:,i));

        for j = 1:nd

            % azimuth and zenith of each body on each sensor frame
            % n.b. boresight direction is along +Z_URF
            q_rbf2urf(:,j) = cp_quat_conj(cp_dcm_to_quat(this(j).dcm0_v));
            q_eci2urf(:,:,j) = cp_quat_mult(q_eci2rbf, q_rbf2urf(:,j));
            sc2body_dir_urf(:,:,i,j) = cp_rotframe(sc2body_dir_eci(:,:,i), q_eci2urf(:,:,j));
            sc2body_sph_urf(:,:,i,j) = cp_sph_coord(sc2body_dir_urf(:,:,i,j));
            sc2body_az_urf(i,:,j) = cp_wrap_to_2pi(sc2body_sph_urf(2,:,i,j));
            sc2body_ze_urf(i,:,j) = pi/2-sc2body_sph_urf(3,:,i,j);

            % fov
            switch p.fov_shape
                case 'circle'
                    fov_ang = this(j).fov_v;
                case 'square'
                    % ...
                    error('To be implemented')
            end

            % baffle
            switch p.baffle_shape
                case 'circle'
                    switch bodies_list{i}
                        case 'earth'
                            baffle_ang = this(j).baffle_earth_v;
                        case 'sun'
                            baffle_ang = this(j).baffle_sun_v;
                        otherwise
                            baffle_ang = nan;
                    end
                case 'square'
                    error('To be implemented')
            end       

            % separation angles required
            sepAngReq_fov(i,:,j) = fov_ang + angDim_body(1,:,i);
            sepAngReq_baffle(i,:,j) =  baffle_ang + angDim_body(1,:,i);

            % blindings
            ix_blind_fov(i,:,j) = sc2body_ze_urf(i,:,j) < sepAngReq_fov(i,:,j);
            ix_blind_baffle(i,:,j) = sc2body_ze_urf(i,:,j) < sepAngReq_baffle(i,:,j);

        end
    end
end
