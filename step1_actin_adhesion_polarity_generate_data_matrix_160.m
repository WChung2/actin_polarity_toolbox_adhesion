%% Generate data matrix

% Load Relion 3D refinement / job006
% [~,particle_indx,angle_rot,angle_tilt,angle_psi,origin_x,origin_y,~,rlnLogLikeliContribution,rlnMaxValueProbDistribution,rlnNrOfSignificantSamples] = actin_polarity_parse_data_file_refine3d('/home/Medalia/Projects7/WChung/actin_adhesion/relion/Refine3D/job082/run_data_mod.star');
[~,particle_indx,angle_rot,angle_tilt,angle_psi,origin_x,origin_y,~,rlnLogLikeliContribution,rlnMaxValueProbDistribution,rlnNrOfSignificantSamples] = actin_polarity_parse_data_file_refine3d_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Refine3D/job019/run_data_mod.star');
% [~,particle_indx,angle_rot,angle_tilt,angle_psi,origin_x,origin_y,~,rlnLogLikeliContribution,rlnMaxValueProbDistribution,rlnNrOfSignificantSamples] = actin_polarity_parse_data_file_refine3d_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Refine3D/job105/run_data_mod.star');

% Sort Relion transformations
[particle_indx_sort,indx_sort] = sort(particle_indx);
angle_rot_sort = angle_rot(1,indx_sort);
angle_tilt_sort = angle_tilt(1,indx_sort);
angle_psi_sort = angle_psi(1,indx_sort);
tx_sort = origin_x(1,indx_sort);
ty_sort = origin_y(1,indx_sort);
rlnLogLikeliContribution_sort = rlnLogLikeliContribution(1,indx_sort);
rlnMaxValueProbDistribution_sort = rlnMaxValueProbDistribution(1,indx_sort);
rlnNrOfSignificantSamples_sort = rlnNrOfSignificantSamples(1,indx_sort);

% Parse data file / Refine3D / job006 / template up
% [~,particle_indx_up,~,~,angle_psi_up] = actin_polarity_parse_data_file_refine3d('/home/Medalia/Projects7/WChung/actin_adhesion/relion/Refine3D/job082/run_data_mod.star');
[~,particle_indx_up,~,~,angle_psi_up] = actin_polarity_parse_data_file_refine3d_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Refine3D/job019/run_data_mod.star');
% [~,particle_indx_up,~,~,angle_psi_up] = actin_polarity_parse_data_file_refine3d_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Refine3D/job105/run_data_mod.star');

% Parse data file / Refine3D / job007 / template down
% [~,particle_indx_down,~,~,angle_psi_down] = actin_polarity_parse_data_file_refine3d('/home/Medalia/Projects7/WChung/actin_adhesion/relion/Refine3D/job086/run_data_mod.star');
[~,particle_indx_down,~,~,angle_psi_down] = actin_polarity_parse_data_file_refine3d_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Refine3D/job023/run_data_mod.star');
% [~,particle_indx_down,~,~,angle_psi_down] = actin_polarity_parse_data_file_refine3d_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Refine3D/job123/run_data_mod.star');

% Check particle index consistency
tom_dev(particle_indx-particle_indx_up);
tom_dev(particle_indx-particle_indx_down);
tom_dev(particle_indx_up-particle_indx_down);

% Sort psi up and down
angle_psi_up_sort = angle_psi_up(1,indx_sort);
angle_psi_down_sort = angle_psi_down(1,indx_sort);

% Calculate segment orientational sensitivity
segment_orientational_sensitivity_pre = abs(angle_psi_up_sort-angle_psi_down_sort);
segment_orientational_sensitivity = (-1).*ones(1,size(segment_orientational_sensitivity_pre,2));
segment_orientational_sensitivity(   find(segment_orientational_sensitivity_pre >= 150 & segment_orientational_sensitivity_pre <= 210)   ) = +1;

% Extract prealignment of good particles
alg_prealg.indx = particle_indx_sort;
alg_prealg.ccc = zeros(1,size(particle_indx_sort,2));
alg_prealg.phi = zeros(1,size(particle_indx_sort,2));
alg_prealg.psi = zeros(1,size(particle_indx_sort,2));
alg_prealg.theta = zeros(1,size(particle_indx_sort,2));
alg_prealg.tx = zeros(1,size(particle_indx_sort,2));
alg_prealg.ty = zeros(1,size(particle_indx_sort,2));
alg_prealg.tz = zeros(1,size(particle_indx_sort,2));

% Add 90 deg rotation
alg_rot.indx = particle_indx_sort;
alg_rot.ccc = zeros(1,size(particle_indx_sort,2));
alg_rot.phi = zeros(1,size(particle_indx_sort,2));
alg_rot.psi = zeros(1,size(particle_indx_sort,2));
alg_rot.theta = zeros(1,size(particle_indx_sort,2)); 
alg_rot.tx = zeros(1,size(particle_indx_sort,2));
alg_rot.ty = zeros(1,size(particle_indx_sort,2));
alg_rot.tz = zeros(1,size(particle_indx_sort,2));

% Combine transformations
alg_com.indx = particle_indx_sort;
alg_com.ccc = zeros(1,size(particle_indx_sort,2));
alg_com.phi = zeros(1,size(particle_indx_sort,2));
alg_com.psi = zeros(1,size(particle_indx_sort,2));
alg_com.theta = zeros(1,size(particle_indx_sort,2));
alg_com.tx = zeros(1,size(particle_indx_sort,2));
alg_com.ty = zeros(1,size(particle_indx_sort,2));
alg_com.tz = zeros(1,size(particle_indx_sort,2));

for k=1:size(alg_com.indx,2)
    
    rotations_tom = [     alg_prealg.phi(k)      alg_prealg.psi(k)      alg_prealg.theta(k);...
                             alg_rot.phi(k)         alg_rot.psi(k)         alg_rot.theta(k)];
             
    translations_tom = [  alg_prealg.tx(k)      alg_prealg.ty(k)      alg_prealg.tz(k);...
                             alg_rot.tx(k)         alg_rot.ty(k)         alg_rot.tz(k)];
    
    [r_sum,t_sum] = actin_polarity_combine_transformations(rotations_tom,translations_tom,angle_psi_sort(k),angle_tilt_sort(k),angle_rot_sort(k),tx_sort(k),ty_sort(k),0);
    
    alg_com.phi(k) = r_sum(1);
    alg_com.psi(k) = r_sum(2);
    alg_com.theta(k) = r_sum(3);
    alg_com.tx(k) = t_sum(1);
    alg_com.ty(k) = t_sum(2);
    alg_com.tz(k) = t_sum(3);
    
end

% Invert combined alignment
alg_com_inv.indx = particle_indx_sort;
alg_com_inv.ccc = zeros(1,size(particle_indx_sort,2));
alg_com_inv.phi = zeros(1,size(particle_indx_sort,2));
alg_com_inv.psi = zeros(1,size(particle_indx_sort,2));
alg_com_inv.theta = zeros(1,size(particle_indx_sort,2));
alg_com_inv.tx = zeros(1,size(particle_indx_sort,2));
alg_com_inv.ty = zeros(1,size(particle_indx_sort,2));
alg_com_inv.tz = zeros(1,size(particle_indx_sort,2));

for k=1:size(alg_com.indx,2)
    
    rotations = [alg_com.phi(k) alg_com.psi(k) alg_com.theta(k)];
             
    translations = [alg_com.tx(k) alg_com.ty(k) alg_com.tz(k)];
    
    [r_sum_inv,t_sum_inv] = InvertTransformations(rotations,translations);
    
    alg_com_inv.phi(k) = r_sum_inv(1);
    alg_com_inv.psi(k) = r_sum_inv(2);
    alg_com_inv.theta(k) = r_sum_inv(3);
    alg_com_inv.tx(k) = t_sum_inv(1);
    alg_com_inv.ty(k) = t_sum_inv(2);
    alg_com_inv.tz(k) = t_sum_inv(3);
    
end

% Construct base path
actin_adhesion_base_path_cells_lessbent;

% Initialize actin data matrix
actin_data_matrix = zeros(2341836,38);

% Loop over all actin particles
zaehler = 1;
zaehler_particles = 1; 

for k=1:size(BasePath,2)
    
    
    % Perform task
    plist = importdata([BasePath{k} 'f_list_lessbent_clear.txt']);
    
    for p=1:size(plist,1) 
        % Fill actin data matrix
        actin_data_matrix(zaehler,1)  = k;                                                      % particle_t
        actin_data_matrix(zaehler,2)  = p;                                                      % particle_p
        actin_data_matrix(zaehler,3)  = zaehler;                                                % particle_n
        
        actin_data_matrix(zaehler,5)  = 4.*plist(p,2);                                          % actin x - coordinate / bin 0 / actin segment position
        actin_data_matrix(zaehler,6)  = 4.*plist(p,3);                                          % actin y - coordinate / bin 0 / actin segment position
        actin_data_matrix(zaehler,7) =  4.*plist(p,4);                                          % actin z - coordinate / bin 0 / actin segment position
        
        actin_data_matrix(zaehler,8) =  4.*plist(p,2);                                          % actin x - coordinate / bin 0 / actin segment position / alignment corrected
        actin_data_matrix(zaehler,9) =  4.*plist(p,3);                                          % actin y - coordinate / bin 0 / actin segment position / alignment corrected
        actin_data_matrix(zaehler,10) = 4.*plist(p,4);                                          % actin z - coordinate / bin 0 / actin segment position / alignment corrected
        
        % Segment to filament assignment
        actin_data_matrix(zaehler,38) = plist(p,1);                                             % filament_index / Segment k-p-zaehler belongs to plist_filaments(p,1)
        
        if ~isempty(find(particle_indx_sort==zaehler))
            
            % Particle reconstructed yes
            actin_data_matrix(zaehler,4)  = 1;                                                  % reconstruced 0/1
            
            
            
            % Measure the correlation between the average of all the particles with an individual particle
            actin_data_matrix(zaehler,35) = 1;
            
           
            % Add segment orientational sensitivity
            actin_data_matrix(zaehler,36) = segment_orientational_sensitivity_pre(1,zaehler_particles);
            actin_data_matrix(zaehler,37) = segment_orientational_sensitivity(1,zaehler_particles);
            
            % Fill actin data matrix
            
            actin_data_matrix(zaehler,8) =  4.*plist(p,2)  +        alg_com.tx(1,zaehler_particles);                              % actin x - coordinate / bin 0 / actin segment position / alignment corrected
            actin_data_matrix(zaehler,9) =  4.*plist(p,3)  +        alg_com.ty(1,zaehler_particles);                              % actin y - coordinate / bin 0 / actin segment position / alignment corrected
            actin_data_matrix(zaehler,10) = 4.*plist(p,4)  +        alg_com.tz(1,zaehler_particles);                              % actin z - coordinate / bin 0 / actin segment position / alignment corrected
            
            actin_data_matrix(zaehler,11) = alg_com.phi(1,zaehler_particles);                                                                   % all alignments / phi
            actin_data_matrix(zaehler,12) = alg_com.psi(1,zaehler_particles);                                                                   % all alignments / psi
            actin_data_matrix(zaehler,13) = alg_com.theta(1,zaehler_particles);                                                                 % all alignments / theta
            
            actin_data_matrix(zaehler,14) = alg_com.tx(1,zaehler_particles);                                                                    % all alignments / tx
            actin_data_matrix(zaehler,15) = alg_com.ty(1,zaehler_particles);                                                                    % all alignments / ty
            actin_data_matrix(zaehler,16) = alg_com.tz(1,zaehler_particles);                                                                    % all alignments / tz
            
            actin_data_matrix(zaehler,17) = sqrt((alg_com.tx(1,zaehler_particles)).^2 + (alg_com.ty(1,zaehler_particles)).^2 + (alg_com.tz(1,zaehler_particles)).^2);    % beta
            
            actin_data_matrix(zaehler,18) = alg_com_inv.phi(1,zaehler_particles);                                                               % all alignments inverted / phi
            actin_data_matrix(zaehler,19) = alg_com_inv.psi(1,zaehler_particles);                                                               % all alignments inverted / psi
            actin_data_matrix(zaehler,20) = alg_com_inv.theta(1,zaehler_particles);                                                             % all alignments inverted / theta
            
            actin_data_matrix(zaehler,21) = alg_com_inv.tx(1,zaehler_particles);                                                                % all alignments inverted / tx
            actin_data_matrix(zaehler,22) = alg_com_inv.ty(1,zaehler_particles);                                                                % all alignments inverted / ty
            actin_data_matrix(zaehler,23) = alg_com_inv.tz(1,zaehler_particles);                                                                % all alignments inverted / tz
            
            % Calculate inverse alignment vector
            [vector_z_rot] = actin_polarity_prepare_z_vector(alg_com_inv.phi(1,zaehler_particles),alg_com_inv.psi(1,zaehler_particles),alg_com_inv.theta(1,zaehler_particles));
            
            actin_data_matrix(zaehler,24) = vector_z_rot(1);                                              % alignment vector / x component
            actin_data_matrix(zaehler,25) = vector_z_rot(2);                                              % alignment vector / y component
            actin_data_matrix(zaehler,26) = vector_z_rot(3);                                              % alignment vector / z component
            
            % Parameter Relion 3D refinement
            
            actin_data_matrix(zaehler,27) = angle_rot_sort(1,zaehler_particles);                          % angle rot
            actin_data_matrix(zaehler,28) = angle_tilt_sort(1,zaehler_particles);                         % angle tilt
            actin_data_matrix(zaehler,29) = angle_psi_sort(1,zaehler_particles);                          % angle psi
            
            actin_data_matrix(zaehler,30) = tx_sort(1,zaehler_particles);                                 % tx
            actin_data_matrix(zaehler,31) = ty_sort(1,zaehler_particles);                                 % ty
            
            actin_data_matrix(zaehler,32) = rlnLogLikeliContribution_sort(1,zaehler_particles);           % rlnLogLikeliContribution
            actin_data_matrix(zaehler,33) = rlnMaxValueProbDistribution_sort(1,zaehler_particles);        % rlnMaxValueProbDistribution
            actin_data_matrix(zaehler,34) = rlnNrOfSignificantSamples_sort(1,zaehler_particles);          % rlnNrOfSignificantSamples
            
            zaehler_particles = zaehler_particles + 1;
            
        end
        
        zaehler = zaehler + 1;
        
        %disp(zaehler);
        
    end
    
    clear plist;
    
    
    
    disp(k);
    
end


% Save workspace
save('/home/Medalia/Projects7/WChung/actin_adhesion/mapping/actin_polarity_mapping_160_step_1_2.mat');
clear all;

