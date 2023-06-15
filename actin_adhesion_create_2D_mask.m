%% Create 2d mask for 2d classification

% Create volumes of zeros
mask = zeros(160,160);

% Fill in ones, 56 pixel = 12.3 nm
% mask(:,80-28+1:80+28) = 1; % 56 pixels have the value of 1
mask(:,80-30+1:80+30) = 1; % 60 pixels have the value of 1



% Create soft edges
mask(:,50:50) = 0.9;
mask(:,110:110) = 0.9;
mask(:,49:49) = 0.8;
mask(:,111:111) = 0.8;
mask(:,48:48) = 0.7;
mask(:,112:112) = 0.7;
mask(:,47:47) = 0.6;
mask(:,113:113) = 0.6;
mask(:,46:46) = 0.5;
mask(:,114:114) = 0.5;
mask(:,45:45) = 0.4;
mask(:,115:115) = 0.4;
mask(:,44:44) = 0.3;
mask(:,116:116) = 0.3;
mask(:,43:43) = 0.2;
mask(:,117:117) = 0.2;
mask(:,42:42) = 0.1;
mask(:,115:118) = 0.1;

%Write to file
tom_mrcwrite(single(mask),'name',['/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/rsc/mask_160_tropo.mrc']);

% % Go to directory of classes
% cd /mnt/ome/Projects7/Actin_MEFs/compare_z_proj_dim/relion_32_newCTF_no_mask/Class2D/job009;
%
% % Read mrcs
% class = tom_mrcread('run_it050_classes.mrcs');
% class = class.Value;
%
% % Extract class #8
% class = class(:,:,8);
%
% % Overlap mask & class for inspection
% merge = mask.*class;