%% bin tomogram to getter better contrast for define substrate

% bin3 of the bin2 tomogram
actin_adhesion_base_path_cells

for i = 1:size(BasePath,2)
    
    % convert mod file to txt file
    cd(BasePath{i})
    tomoname = dir(['*','.rawtlt']);
    [~,name,~] = fileparts(tomoname.name);
    unix(['binvol -i ' name '_rec.mrc -o tomo_bin.mrc -x 3 -y 3 -z 3']);
    clearvars name tomoname
    
end



%%  extract base and the normal angle

actin_adhesion_base_path_cells_lessbent;


for i = 1:size(BasePath,2)
    
    % convert mod file to txt file
    cd(BasePath{i})
    unix(['model2point b_list b_list.txt']);
    
end



% get normal angle of the base
actin_adhesion_base_path_cells_lessbent;

for i = 1:size(BasePath,2)
    cd(BasePath{i})
    bin_bPlist = importdata('b_list.txt');
    bPlist = bin_bPlist.*3;
    corbPlist = bPlist;
    [n,V,p] = affine_fit(corbPlist);
    Zv = [0 0 1];
    n = n';
    corR = vrrotvec(n,Zv);
    rotm = axang2rotm(corR);
    eul = rotm2eul(rotm,'XYZ');
    corb_angle(i,:) = eul;
    clear bPlist;
    disp(i)

    
end

clearvars -except corb_angle


save('/home/Medalia/Projects7/WChung/actin_adhesion/mapping/corb_angle.mat');
clear

%% extract adhesion vector

% get direction from 3dmod, 4st to 3rd, and 2nd to 1st cor from v_list
% if v_wrong exist, the direction is not defined


actin_adhesion_base_path_cells_lessbent;

for k = 1:size(BasePath,2)
    
    cd(BasePath{k})
    
    if isfile('v_list')
        unix(['model2point v_list v_list.txt']);
        
    elseif isfile('v_wrong')
        
        
    else
        disp([num2str(k) " something is wrong"])
        
    end
    
    
end

corv_angle = {};

for k = 1:size(BasePath,2)
    cd(BasePath{k})
    
    if isfile('v_list.txt')
        
        
        
        VPlist = importdata('v_list.txt');
        if size(VPlist,1) == 2
            corv_angle{k,1} = VPlist(2,:)-VPlist(1,:);
            
        elseif size(VPlist,1) == 4
            corv_angle{k,1} = VPlist(2,:)-VPlist(1,:);
            corv_angle{k,2} = VPlist(4,:)-VPlist(3,:);
            
        else
            disp([num2str(k) " vlist is wrong"])
            
            
            
        end
        
        
        clear VPlist
        disp(k)
        
    elseif isfile('v_wrong')
        
    else
        
        disp([num2str(k) " something is wrong"])
        
        
    end
    
    
end

save('/home/Medalia/Projects7/WChung/actin_adhesion/mapping/corv_angle.mat');

clear