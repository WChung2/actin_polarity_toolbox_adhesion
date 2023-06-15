% Actin polarity / averaging / box size 104 pixels 

%% Create particle list and alignment structure

% Generate Particle BasePath
% actin_adhesion_particle_base_path_cells;
% actin_adhesion_base_path_cells;

actin_adhesion_particle_base_path_cells_lessbent;
actin_adhesion_base_path_cells_lessbent;


for k=1:size(P_BasePath,2)
    
    %read f_list
    
    f_list = importdata([BasePath{k} 'f_list_lessbent_clear.txt']);
    
    
    particlename = dir([P_BasePath{k} 'imod_subtomo_extract/']);
    
   
        
    parfor p=1:size(f_list,1)
        
        % Identify particle
        actin_particle_pre = [P_BasePath{k} 'imod_subtomo_extract/' particlename(p+2).name];
        disp(actin_particle_pre);
        
        % Create new particle identifier
        actin_particle = [P_BasePath{k} 'actin_particle_' actin_polarity_construct_particle_identifier_adhesion(k,f_list(p,1),p) '.mrc'];
        disp(actin_particle);
        
        % Rename particle
        unix(['mv ' actin_particle_pre ' ' actin_particle]);
        
        
    end
         
    
    clear f_list;
    
    
end

clear

% Build particle list
% actin_adhesion_particle_base_path_cells;
% actin_adhesion_base_path_cells;
actin_adhesion_particle_base_path_cells_lessbent;
actin_adhesion_base_path_cells_lessbent;



zaehler = 1;
particles_list = {};

for k = 1:size(P_BasePath,2)
    
    particlename = dir([P_BasePath{k}]);
    
    for i = 3:size(particlename,1)-1
    
        particles_list{zaehler} = [P_BasePath{k} particlename(i).name];
        
        zaehler = zaehler + 1;
    
    end
    disp(k)

end




zaehler = 1;
for k=1:size(particles_list,2)
    
    particles_list_combined{zaehler} = particles_list{k};
    zaehler = zaehler + 1;
    
    particles_list_combined{zaehler} = 'singleaxiswedge 30 80';
    zaehler = zaehler + 1;
    
    particles_list_combined{zaehler} = 'allpass';
    zaehler = zaehler + 1;
    
    particles_list_combined{zaehler} = ' ';
    zaehler = zaehler + 1;
    
end

% Write empty file
placeh = num2str(0);
save('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/actin_polarity_particles_list_bin0_160_lessbent.txt', 'placeh', '-ASCII');

% Write particles into file
fid = fopen('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/actin_polarity_particles_list_bin0_160_lessbent.txt', 'w');
for k=1:size(particles_list_combined,2)
    fprintf(fid, '%s\n', char(particles_list_combined{k}));        
end
fclose(fid);

% Save workspace
save('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/workspace_actin_polarity_averaging_160_step_1.mat');
clear all;


%% Project particles

% Parse particles list
[particles_list] = ParseParticleList('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/actin_polarity_particles_list_bin0_160_lessbent.txt');

% Make projection folder

% projectPath = '/home/Medalia/Projects5/Wen-Lu/actin_adhesion_160/Particles_Proj_160_z50/';
projectPath = '/home/Medalia/Projects7/WChung/actin_adhesion/Particles_Proj_160_z50/';
mkdir(projectPath)

% cd /home/Medalia/Projects7/WChung/actin_MEF_tomogram/Manual_Segmentation/20190520_VV/Particles_Proj_160_z50_clear/



% Project particles
parfor k=1:size(particles_list,2)
    
    % Load particle
    particle = double(tom_mrcreadf(particles_list{k}));
    
    % Project actin filament
    proj = mean(particle(:,:,80-25+1:80+25),3);%---> Projection thickness = 11 nm (50 pixels)
    
    % Normalize particle
    proj = (proj - mean(proj(:)))./std(proj(:));
    
    % Dissect file name into parts
    [~,particle_name] = fileparts(particles_list{k});
    
    % Write to disk
    tom_mrcwrite(single(proj),'name',[projectPath particle_name '.mrc']);
    
                      % Prepare proj list
                      proj_list{k} = [projectPath particle_name '.mrc' '   ' projectPath  particle_name(1:21) '.mrc'];
                      
    disp([num2str(k/size(particles_list,2)*100) ' %']);
    
end

% Write proj list
placeh = num2str(0);
       save('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/actin_polarity_particles_list_bin0_160_proj_lessbent.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/actin_polarity_particles_list_bin0_160_proj_lessbent.star', 'w');

% Add the following header:
% data_
% loop_
% _rlnImageName
% _rlnMicrographName
str1=["data_"];
str2=["loop_"];
str3=["_rlnImageName #1"];
str4=["_rlnMicrographName #2"];

fprintf(fid, '%s\n', str1, str2, str3, str4);
for k=1:size(proj_list,2)
    fprintf(fid, '%s\n', char(proj_list{k}));
end
fclose(fid);



% Normalize with Relion preprocess
unix(['mpirun -np 64 relion_preprocess_mpi --operate_on /home/Medalia/Projects7/WChung/actin_adhesion/averaging/actin_polarity_particles_list_bin0_160_proj_lessbent.star --norm --bg_radius 80 --invert_contrast --operate_out /home/Medalia/Projects7/WChung/actin_adhesion/averaging/actin_polarity_particles_list_bin0_160_proj_lessbent_norm_inv --set_angpix 2.20651']);

% Save workspace
save('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/workspace_actin_polarity_averaging_160_step_2.mat');
% clear all;



%% Prepare 2D classification with priors / Class2D / job001 / it010 / prealignment with template library

% Parse data file / "this is the prealignment"
[alg_struct] = actin_polarity_parse_data_file_class2d_prealg_adhesion('/home/Medalia/Projects7/WChung/actin_adhesion/relion/Class2D/job193/run_it009_data_mod.star');


pixelsize = 2.20651;

% Create new starfile
relion_data_file_out = {};
zaehler = 1;
for k=1:size(alg_struct,2)
    
    if ~ismember(alg_struct(k).group_number,[29  111  112])
        relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' alg_struct(k).micrograph_name '   ' num2str(0) '   ' num2str(0) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(1) '   ' num2str((alg_struct(k).origin_x)*pixelsize) '   ' num2str((alg_struct(k).origin_y)*pixelsize)];
        
        zaehler = zaehler + 1;
    end
end

% Write starfile
placeh = num2str(0);
       save('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/particles_job193_it009_prior_edit.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/particles_job193_it009_prior_edit.star', 'w');
% Add the following header:
% data_
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnOriginXPrior #6 
% _rlnOriginYPrior #7 

strB = [""];
str1=["# version 30001"];
str2=["data_optics"];
str3=["loop_"];
str4=["_rlnOpticsGroup #1"];
str5=["_rlnOpticsGroupName #2"];
str6=["_rlnImagePixelSize #3"];
str7=["_rlnImageSize #4"];
str34=["_rlnImageDimensionality #5"];
str8=["_rlnVoltage #6"];
str9=["_rlnSphericalAberration #7"];
str10=["_rlnAmplitudeContrast #8"];
str11=["           1 opticsGroup1   2.20651    160       2   300  2.7   0.07"];
str12=["# version 30001"];
str13=["data_particles"];
str14=["loop_"];
str15=["_rlnImageName #1"];
str16=["_rlnMicrographName #2"];
str17=["_rlnAngleRotPrior #3"];
str18=["_rlnAngleTiltPrior #4"];
str19=["_rlnAnglePsiPrior #5"];
str20=["_rlnOpticsGroup #6"];
str21=["_rlnOriginXPriorAngst #7"];
str22=["_rlnOriginYPriorAngst #8"];

fprintf(fid, '%s\n', str1, strB, str2, strB, str3, str4, str5, str6, str7, str34, str8, str9, str10, str11, strB, strB, str12, strB, str13, strB, str14, str15, str16, str17, str18, str19, str20, str21, str22);
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);

% Add header manually

% data_
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnOriginXPrior #6 
% _rlnOriginYPrior #7 


strB = [""];
str1=["# version 30001"];
str2=["data_optics"];
str3=["loop_"];
str4=["_rlnOpticsGroup #1"];
str5=["_rlnOpticsGroupName #2"];
str6=["_rlnImagePixelSize #3"];
str7=["_rlnImageSize #4"];
str34=["_rlnImageDimensionality #5"];
str8=["_rlnVoltage #6"];
str9=["_rlnSphericalAberration #7"];
str10=["_rlnAmplitudeContrast #8"];
str11=["           1 opticsGroup1   2.2065    160       2   300  2.7   0.07"];
str12=["# version 30001"];
str13=["data_particles"];
str14=["loop_"];
str15=["_rlnImageName #1"];
str16=["_rlnMicrographName #2"];
str17=["_rlnAngleRotPrior #3"];
str18=["_rlnAngleTiltPrior #4"];
str19=["_rlnAnglePsiPrior #5"];
str20=["_rlnGroupNumber #6"];
str21=["_rlnAngleRot #7"];
str22=["_rlnAngleTilt #8"];
str23=["_rlnAnglePsi #9"];
str24=["_rlnClassNumber #10"];
str25=["_rlnNormCorrection #11"];
str26=["_rlnLogLikeliContribution #12"];
str27=["_rlnMaxValueProbDistribution #13"];
str28=["_rlnNrOfSignificantSamples #14"];
str29=["_rlnOpticsGroup #15"];
str30=["_rlnOriginXAngst #16"];
str31=["_rlnOriginYAngst #17"];
str32=["_rlnOriginXPriorAngst #18"];
str33=["_rlnOriginYPriorAngst #19"];


fprintf(fid, '%s\n', str1, strB, str2, strB, str3, str4, str5, str6, str7, str34, str8, str9, str10, str11, strB, strB, str12, strB, str13, strB, str14, str15, str16, str17, str18, str19, str20, str21, str22, str23, str24, str25, str26, str27, str28, str29, str30, str31, str32, str33);
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);

%% prepare star file for relion4 2d classification


pixelsize = 2.2065;

[alg_struct] = actin_polarity_parse_data_file_class2d_finealg_adhesion('/home/Medalia/Projects7/WChung/actin_adhesion/relion/Select/job196/particles_mod.star');

% Create new starfile
zaehler = 1;
for k=1:size(alg_struct,2)
    
    relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' alg_struct(k).micrograph_name '   '  num2str(0) '   ' num2str(0) '   ' num2str(alg_struct(k).AnglePsiPrior) '   ' num2str(alg_struct(k).group_number) '   ' ...
       num2str(0) '   ' num2str(0) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(alg_struct(k).class_indx) '   ' num2str(alg_struct(k).norm_corr) '   ' num2str(alg_struct(k).log_likeli_con) '   ' num2str(alg_struct(k).max_val_prob_dis) '   ' ...
       num2str(alg_struct(k).nr_sig_samples) '   ' num2str(1) '   ' num2str((alg_struct(k).OriginXPrior)*pixelsize) '   ' num2str((alg_struct(k).OriginYPrior)*pixelsize) '   ' num2str((alg_struct(k).origin_x)*pixelsize) '   ' num2str((alg_struct(k).origin_y)*pixelsize)];
    
    zaehler = zaehler + 1;
    
end

% Write starfile
placeh = num2str(0);
       save('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/particles_relion4_select_job196_it009.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/particles_relion4_select_job196_it009.star', 'w');



% Add the following header:
% data_
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnOriginXPrior #6 
% _rlnOriginYPrior #7 

str1=["data_"];
str2=["loop_"];
str3=["_rlnImageName #1"];
str4=["_rlnMicrographName #2"];
str5=["_rlnAngleRotPrior #3"];
str6=["_rlnAngleTiltPrior #4"];
str7=["_rlnAnglePsiPrior #5"];
str8=["_rlnOriginXPrior #6"];
str9=["_rlnOriginYPrior #7"];

fprintf(fid, '%s\n', str1, str2, str3, str4, str5, str6, str7, str8, str9);
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);




% Add header manually



% # version 30001
% 
% data_optics
% 
% loop_ 
% _rlnOpticsGroup #1 
% _rlnOpticsGroupName #2 
% _rlnImagePixelSize #3
% _rlnImageSize #4 
% _rlnImageDimensionality #5 
% _rlnVoltage #6
% _rlnSphericalAberration #7
% _rlnAmplitudeContrast #8 
%            1 opticsGroup1   2.2065    160       2   300  2.7   0.07
%  
% 
% # version 30001
% 
% data_particles
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnGroupNumber #6 
% _rlnAngleRot #7 
% _rlnAngleTilt #8 
% _rlnAnglePsi #9 
% _rlnClassNumber #10 
% _rlnNormCorrection #11 
% _rlnLogLikeliContribution #12 
% _rlnMaxValueProbDistribution #13 
% _rlnNrOfSignificantSamples #14 
% _rlnOpticsGroup #15 
% _rlnOriginXAngst #16 
% _rlnOriginYAngst #17 
% _rlnOriginXPriorAngst #18 
% _rlnOriginYPriorAngst #19 


%% Prepare 3D reconstruction with priors / Class2D / job003 / it100 / Select / job004 / alignment with prior psi angle

% Parse data file / "this is the psi finealignment"
% [alg_struct] = actin_polarity_parse_data_file_class2d_finealg_adhesion('/home/Medalia/Projects7/WChung/actin_adhesion/relion/Select/job126/particles_mod.star');
% [alg_struct] = actin_polarity_parse_data_file_class2d_finealg_adhesion('/home/Medalia/Projects7/WChung/actin_adhesion/relion/Select/job147/particles_mod.star');
% [alg_struct] = actin_polarity_parse_data_file_class2d_finealg_adhesion_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Select/job005/particles_mod.star');
[alg_struct] = actin_polarity_parse_data_file_class2d_finealg_adhesion_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Select/job101/particles_mod.star');
[particles_list] = ParseParticleList('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/actin_polarity_particles_list_bin0_160_lessbent.txt');


particleN = [];
tomorawN = [];
parfor k=1:size(particles_list,2)
    [~,particle_name] = fileparts(particles_list{k});
    name = particle_name;
    idx1 = regexp(name, '_f_', 'end');
    idx2 = regexp(name, '_p_', 'start');
    idx3 = regexp(name, '_t_', 'end');
    idx4 = regexp(name, '_f_', 'start');
    particleN(k,1)= str2double(name(idx1+1 : idx2-1));%filament number
    tomorawN(k,1)= str2double(name(idx3+1 : idx4-1));%tomogram number
    disp(k)
end




zaehler = 1;
counter1 = 0;
counter2 = 0;
for k=1:size(alg_struct,2)
    filamentN = particleN(alg_struct(k).particle_indx,1);
    tomoN = tomorawN(alg_struct(k).particle_indx,1);
    
    if rem(tomoN,2)==1
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(tomoN) '   ' num2str(0) '   ' num2str(90) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(1) '   ' num2str(1)];
         counter1 = counter1 + 1;
    
    elseif rem(tomoN,2)==0
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(tomoN) '   ' num2str(0) '   ' num2str(90) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(2) '   ' num2str(1)];
         counter2 = counter2 + 1;
    end
    
    zaehler = zaehler + 1;
    disp(k)
end





% # version 30001
% 
% data_optics
% 
% loop_ 
% _rlnOpticsGroup #1 
% _rlnOpticsGroupName #2 
% _rlnImagePixelSize #3 
% _rlnImageSize #4 
% _rlnImageDimensionality #5 
% _rlnVoltage #6 
% _rlnSphericalAberration #7 
% _rlnAmplitudeContrast #8 
%            1 opticsGroup1     2.206500          160            2   300.000000     2.700000     0.070000 
%  
% 
% # version 30001
% 
% data_particles
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnGroupNumber #6 
% _rlnAngleRot #7 
% _rlnAngleTilt #8 
% _rlnAnglePsi #9 
% _rlnClassNumber #10 
% _rlnNormCorrection #11 
% _rlnLogLikeliContribution #12 
% _rlnMaxValueProbDistribution #13 
% _rlnNrOfSignificantSamples #14 
% _rlnOpticsGroup #15 
% _rlnOriginXAngst #16 
% _rlnOriginYAngst #17 
% _rlnOriginXPriorAngst #18 
% _rlnOriginYPriorAngst #19 

placeh = num2str(0);
       save('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/rsc/particles_select_relion4_job101.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/rsc/particles_select_relion4_job101.star', 'w');

strB = [""];
str1=["# version 30001"];
str2=["data_optics"];
str3=["loop_"];
str4=["_rlnOpticsGroup #1"];
str5=["_rlnOpticsGroupName #2"];
str6=["_rlnImagePixelSize #3"];
str7=["_rlnImageSize #4"];
str8=["_rlnImageDimensionality #5"];
str9=["_rlnVoltage #6"];
str10=["_rlnSphericalAberration #7"];
str11=["_rlnAmplitudeContrast #8"];
str12=["           1 opticsGroup1   2.2065    160       2   300  2.7   0.07"];
str13=["# version 30001"];
str14=["data_particles"];
str15=["loop_"];
str16=["_rlnImageName #1"];
str17=["_rlnMicrographName #2"];
str18=["_rlnAngleRotPrior #3"];
str19=["_rlnAngleTiltPrior #4"];
str20=["_rlnAnglePsiPrior #5"];
str21=["_rlnRandomSubset #6"];
str22=["_rlnOpticsGroup #7"];


fprintf(fid, '%s\n', str1, strB, str2, strB, str3, str4, str5, str6, str7, str8, str9, str10, str11, str12, strB, strB, str13, strB, str14, strB, str15, str16, str17, str18, str19, str20, str21, str22);
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);







%% prepare tropomyosin selection

% prepare particle list for latter
% [particles_list] = ParseParticleList('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/actin_polarity_particles_list_bin0_160_lessbent.txt');
% particleN = [];
% tomorawN = [];
% parfor k=1:size(particles_list,2)
%     [~,particle_name] = fileparts(particles_list{k});
%     name = particle_name;
%     idx1 = regexp(name, '_f_', 'end');
%     idx2 = regexp(name, '_p_', 'start');
%     idx3 = regexp(name, '_t_', 'end');
%     idx4 = regexp(name, '_f_', 'start');
%     idx5 = regexp(name, '_p_', 'end');
%     particleN(k,1)= str2double(name(idx1+1 : idx2-1));%filament number
%     tomorawN(k,1)= str2double(name(idx3+1 : idx4-1));%tomogram number
%     numbern(k,1) = str2double(name(idx5+1 : end));%particle number
%     disp(k)
% end
% 
% P_n_F_T_list = zeros(size(particleN,1),4);
% 
% parfor k = 1: size(particleN,1)
%     
%     P_n_F_T_list(k,:)= [k numbern(k) particleN(k) tomorawN(k)];
%     
%     
%     
% end
% 
% writematrix(P_n_F_T_list,'/home/Medalia/Projects7/WChung/actin_adhesion/averaging/P_n_F_T_list','Delimiter','tab')






% create star file contained all the selected "filament" 

[alg_struct2] = actin_polarity_parse_data_file_class2d_finealg_adhesion_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Select/job125/particles_mod.star');
P_n_F_T_list = importdata('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/P_n_F_T_list.txt');


zaehler = 1;
numpool = [];
for k = 1: size(alg_struct2,2)
    numpool_temp = [];
    if ~ismember(alg_struct2(k).particle_indx,numpool)
        filamentINDEX = P_n_F_T_list(alg_struct2(k).particle_indx,3);
        tomoINDEX = P_n_F_T_list(alg_struct2(k).particle_indx,4);
        numpool_temp = P_n_F_T_list((filamentINDEX==P_n_F_T_list(:,3) & tomoINDEX==P_n_F_T_list(:,4)),1);



        
        numpool = [numpool;numpool_temp];
        
    
    end
    
    disp(k)
    
end


% Parse data file / "this is the prealignment"
% [alg_struct] = actin_polarity_parse_data_file_class2d_finealg_adhesion('/home/Medalia/Projects7/WChung/actin_adhesion/relion/Select/job196/particles_mod.star');
[alg_struct] = actin_polarity_parse_data_file_class2d_finealg_adhesion_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Select/job101/particles_mod.star');
relion_data_file_out = {};
zaehler = 1;
for k=1:size(alg_struct,2)
    if ismember(alg_struct(k).particle_indx,numpool)
    
        relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' alg_struct(k).micrograph_name '   '  num2str(0) '   ' num2str(0) '   ' num2str(alg_struct(k).AnglePsiPrior) '   ' num2str(alg_struct(k).group_number) '   ' ...
            num2str(0) '   ' num2str(0) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(alg_struct(k).class_indx) '   ' num2str(alg_struct(k).norm_corr) '   ' num2str(alg_struct(k).log_likeli_con) '   ' num2str(alg_struct(k).max_val_prob_dis) '   ' ...
            num2str(alg_struct(k).nr_sig_samples) '   ' num2str(1) '   ' num2str(alg_struct(k).OriginXAngst) '   ' num2str(alg_struct(k).OriginYAngst) '   ' num2str(alg_struct(k).OriginXPriorAngst) '   ' num2str(alg_struct(k).OriginYPriorAngst)];
        
        zaehler = zaehler + 1;
    
    end
    disp(k)
    
end

% Write starfile
placeh = num2str(0);
       save('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/particles_relion4_select_job125_tropo.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/particles_relion4_select_job125_tropo.star', 'w');

strB = [""];
str1=["# version 30001"];
str2=["data_optics"];
str3=["loop_"];
str4=["_rlnOpticsGroup #1"];
str5=["_rlnOpticsGroupName #2"];
str6=["_rlnImagePixelSize #3"];
str7=["_rlnImageSize #4"];
str34=["_rlnImageDimensionality #5"];
str8=["_rlnVoltage #6"];
str9=["_rlnSphericalAberration #7"];
str10=["_rlnAmplitudeContrast #8"];
str11=["           1 opticsGroup1   2.2065    160       2   300  2.7   0.07"];
str12=["# version 30001"];
str13=["data_particles"];
str14=["loop_"];
str15=["_rlnImageName #1"];
str16=["_rlnMicrographName #2"];
str17=["_rlnAngleRotPrior #3"];
str18=["_rlnAngleTiltPrior #4"];
str19=["_rlnAnglePsiPrior #5"];
str20=["_rlnGroupNumber #6"];
str21=["_rlnAngleRot #7"];
str22=["_rlnAngleTilt #8"];
str23=["_rlnAnglePsi #9"];
str24=["_rlnClassNumber #10"];
str25=["_rlnNormCorrection #11"];
str26=["_rlnLogLikeliContribution #12"];
str27=["_rlnMaxValueProbDistribution #13"];
str28=["_rlnNrOfSignificantSamples #14"];
str29=["_rlnOpticsGroup #15"];
str30=["_rlnOriginXAngst #16"];
str31=["_rlnOriginYAngst #17"];
str32=["_rlnOriginXPriorAngst #18"];
str33=["_rlnOriginYPriorAngst #19"];


fprintf(fid, '%s\n', str1, strB, str2, strB, str3, str4, str5, str6, str7, str34, str8, str9, str10, str11, strB, strB, str12, strB, str13, strB, str14, str15, str16, str17, str18, str19, str20, str21, str22, str23, str24, str25, str26, str27, str28, str29, str30, str31, str32, str33);
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);


% Add header manually



% # version 30001
% 
% data_optics
% 
% loop_ 
% _rlnOpticsGroup #1 
% _rlnOpticsGroupName #2 
% _rlnImagePixelSize #3
% _rlnImageSize #4 
% _rlnImageDimensionality #5 
% _rlnVoltage #6
% _rlnSphericalAberration #7
% _rlnAmplitudeContrast #8 
%            1 opticsGroup1   2.2065    160       2   300  2.7   0.07
%  
% 
% # version 30001
% 
% data_particles
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnGroupNumber #6 
% _rlnAngleRot #7 
% _rlnAngleTilt #8 
% _rlnAnglePsi #9 
% _rlnClassNumber #10 
% _rlnNormCorrection #11 
% _rlnLogLikeliContribution #12 
% _rlnMaxValueProbDistribution #13 
% _rlnNrOfSignificantSamples #14 
% _rlnOpticsGroup #15 
% _rlnOriginXAngst #16 
% _rlnOriginYAngst #17 
% _rlnOriginXPriorAngst #18 
% _rlnOriginYPriorAngst #19 

%% Prepare 3D reconstruction with priors for tropo

% Parse data file / "this is the psi finealignment"
% [alg_struct] = actin_polarity_parse_data_file_class2d_finealg_adhesion('/home/Medalia/Projects7/WChung/actin_adhesion/relion/Select/job126/particles_mod.star');
% [alg_struct] = actin_polarity_parse_data_file_class2d_finealg_adhesion('/home/Medalia/Projects7/WChung/actin_adhesion/relion/Select/job147/particles_mod.star');
[alg_struct] = actin_polarity_parse_data_file_class2d_finealg_adhesion_relion4('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Select/job068/particles_mod.star');
[particles_list] = ParseParticleList('/home/Medalia/Projects7/WChung/actin_adhesion/averaging/actin_polarity_particles_list_bin0_160_lessbent.txt');


particleN = [];
tomorawN = [];
parfor k=1:size(particles_list,2)
    [~,particle_name] = fileparts(particles_list{k});
    name = particle_name;
    idx1 = regexp(name, '_f_', 'end');
    idx2 = regexp(name, '_p_', 'start');
    idx3 = regexp(name, '_t_', 'end');
    idx4 = regexp(name, '_f_', 'start');
    particleN(k,1)= str2double(name(idx1+1 : idx2-1));%filament number
    tomorawN(k,1)= str2double(name(idx3+1 : idx4-1));%tomogram number
    disp(k)
end




zaehler = 1;
counter1 = 0;
counter2 = 0;
for k=1:size(alg_struct,2)
    filamentN = particleN(alg_struct(k).particle_indx,1);
    tomoN = tomorawN(alg_struct(k).particle_indx,1);
    
    if rem(tomoN,2)==1
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(tomoN) '   ' num2str(0) '   ' num2str(90) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(1) '   ' num2str(1)];
         counter1 = counter1 + 1;
    
    elseif rem(tomoN,2)==0
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(tomoN) '   ' num2str(0) '   ' num2str(90) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(2) '   ' num2str(1)];
         counter2 = counter2 + 1;
    end
    
    zaehler = zaehler + 1;
    disp(k)
end




% 
% # version 30001
% 
% data_optics
% 
% loop_ 
% _rlnOpticsGroup #1 
% _rlnOpticsGroupName #2 
% _rlnImagePixelSize #3 
% _rlnImageSize #4 
% _rlnImageDimensionality #5 
% _rlnVoltage #6 
% _rlnSphericalAberration #7 
% _rlnAmplitudeContrast #8 
%            1 opticsGroup1     2.206500          160            2   300.000000     2.700000     0.070000 
%  
% 
% # version 30001
% 
% data_particles
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnGroupNumber #6 
% _rlnAngleRot #7 
% _rlnAngleTilt #8 
% _rlnAnglePsi #9 
% _rlnClassNumber #10 
% _rlnNormCorrection #11 
% _rlnLogLikeliContribution #12 
% _rlnMaxValueProbDistribution #13 
% _rlnNrOfSignificantSamples #14 
% _rlnOpticsGroup #15 
% _rlnOriginXAngst #16 
% _rlnOriginYAngst #17 
% _rlnOriginXPriorAngst #18 
% _rlnOriginYPriorAngst #19 

placeh = num2str(0);
       save('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/rsc/particles_select_relion4_job068.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/rsc/particles_select_relion4_job068.star', 'w');
str1=["data_particles"];
str2=["loop_"];
str3=["_rlnImageName #1"];
str4=["_rlnMicrographName #2"];
str5=["_rlnAngleRotPrior #3"];
str6=["_rlnAngleTiltPrior #4"];
str7=["_rlnAnglePsiPrior #5"];
str8=["_rlnRandomSubset #6"];
str9=["_rlnOpticsGroup #7"];


fprintf(fid, '%s\n', str1, str2, str3, str4, str5, str6, str7, str8, str9);
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);





%% Prepare local alignment and classification based on Refine3D / job006

% Parse data file
[alg_struct] = actin_polarity_parse_data_file_class3d('/home/Medalia/Projects7/WChung/actin_MEF_tomogram/Manual_Segmentation/20190520_VV/relion_proj/Refine3D/job282/run_data_mod.star');

% Create new starfile
zaehler = 1;
for k=1:size(alg_struct,2)
        
    relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' alg_struct(k).micrograph_name '   ' num2str(alg_struct(k).angle_rot) '   ' num2str(alg_struct(k).angle_tilt) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(alg_struct(k).RandomSubset)];% '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];

    
    zaehler = zaehler + 1;
%     disp(k)
end



% zaehler = 1;
% group_indx = 1;
% for k=1:size(alg_struct,2)
%     
%     if group_indx == 1
%          relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(group_indx) '   ' num2str(alg_struct(k).angle_rot) '   ' num2str(alg_struct(k).angle_tilt) '   ' num2str(alg_struct(k).angle_psi)];% '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];
%          group_indx = 2;
%     else
%          relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(group_indx) '   ' num2str(alg_struct(k).angle_rot) '   ' num2str(alg_struct(k).angle_tilt) '   ' num2str(alg_struct(k).angle_psi)];% '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];
%          group_indx = 1;
%     end
%     
%     zaehler = zaehler + 1;
%     
% end

% Write starfile and add header manually
placeh = num2str(0);
       save('/home/Medalia/Projects7/WChung/actin_MEF_tomogram/Manual_Segmentation/20190520_VV/relion_proj/rsc/particles_job282_prior_mod.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/WChung/actin_MEF_tomogram/Manual_Segmentation/20190520_VV/relion_proj/rsc/particles_job282_prior_mod.star', 'w');
str1=["data_"];
str2=["loop_"];
str3=["_rlnImageName #1"];
str4=["_rlnMicrographName #2"];
str5=["_rlnAngleRotPrior #3"];
str6=["_rlnAngleTiltPrior #4"];
str7=["_rlnAnglePsiPrior #5"];
str8=["_rlnRandomSubset #6"];

% str8=["_rlnOriginXPrior #6"];
% str9=["_rlnOriginYPrior #7"];


fprintf(fid, '%s\n', str1, str2, str3, str4, str5, str6, str7, str8);%, str8, str9);
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);

% Add header manually

% data_
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnOriginXPrior #6 
% _rlnOriginYPrior #7 

%% Prepare local alignment and classification based on Refine3D / job006

% Parse data file
[alg_struct] = actin_polarity_parse_data_file_class3d_2('/home/Medalia/Projects7/WChung/actin_MEF_tomogram/Manual_Segmentation/20190520_VV/relion_proj/Refine3D/job177/run_data_mod.star');

% Create new starfile
zaehler = 1;
group_indx = 1;
for k=1:size(alg_struct,2)
    
    if group_indx == 1
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(group_indx) '   ' num2str(alg_struct(k).angle_rot) '   ' num2str(alg_struct(k).angle_tilt) '   ' num2str(alg_struct(k).angle_psi)];% '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];
         group_indx = 2;
    else
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(group_indx) '   ' num2str(alg_struct(k).angle_rot) '   ' num2str(alg_struct(k).angle_tilt) '   ' num2str(alg_struct(k).angle_psi)];% '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];
         group_indx = 1;
    end
    
    zaehler = zaehler + 1;
    
end

% Write starfile and add header manually
placeh = num2str(0);
       save('/home/Medalia/Projects7/WChung/actin_MEF_tomogram/Manual_Segmentation/20190520_VV/relion_proj/rsc/particles_job177_prior_mod.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/WChung/actin_MEF_tomogram/Manual_Segmentation/20190520_VV/relion_proj/rsc/particles_job177_prior_mod.star', 'w');
str1=["data_"];
str2=["loop_"];
str3=["_rlnImageName #1"];
str4=["_rlnMicrographName #2"];
str5=["_rlnAngleRotPrior #3"];
str6=["_rlnAngleTiltPrior #4"];
str7=["_rlnAnglePsiPrior #5"];
% str8=["_rlnOriginXPrior #6"];
% str9=["_rlnOriginYPrior #7"];


fprintf(fid, '%s\n', str1, str2, str3, str4, str5, str6, str7);%, str8, str9);
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);

% Add header manually

% data_
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnOriginXPrior #6 
% _rlnOriginYPrior #7 