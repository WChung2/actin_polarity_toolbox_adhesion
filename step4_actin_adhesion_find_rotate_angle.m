% Create 10 random 3d points for explanation
% cor = round(rand(10,3).* 1024);% This are your segment coordinates eg. dots

%% write rotation psi and tilt angle into filament structure
load('/home/Medalia/Projects7/WChung/actin_adhesion/mapping/corv_angle.mat','corv_angle');
load('/home/Medalia/Projects7/WChung/actin_adhesion/mapping/corb_angle.mat','corb_angle');

%Load filaments and BasePath
actin_adhesion_base_path_cells_lessbent;

for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:189
    FilamentStruct2{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step4/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step4/

for k=[1:28 30:110 113:189]
    
    % load filament structure
    load(FilamentStruct{k},'filament_struct_expanded');
    
   
    for i = 1:size(filament_struct_expanded,2)
        filament_struct_expanded(i).tomogram_stressF_vector = [corv_angle{k,1};corv_angle{k,2}];
        filament_struct_expanded(i).tomogram_base_corr =corb_angle(k,:);
                
    end
    save(FilamentStruct2{k},'filament_struct_expanded');
    disp(k) 
    

end
clear