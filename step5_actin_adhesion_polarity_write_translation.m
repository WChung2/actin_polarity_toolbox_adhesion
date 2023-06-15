%% rotate filaments and membrane to psi and base

%Load filaments and BasePath
actin_adhesion_base_path_cells_lessbent;

for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step4/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:189
    FilamentStruct2{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step5/filaments_160_tomo_' num2str(t) '.mat'];
end


mkdir /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step5/


for k = [1:28 30:110 113:189]

    % load filament structure and membrane and base
    load(FilamentStruct{k},'filament_struct_expanded');
    bin_bPlist = importdata([BasePath{k} 'b_list.txt']);
    bPlist = bin_bPlist.*3;
       
    tomo_rot_psi = 0;
    tomo_rot_b = filament_struct_expanded(1).tomogram_base_corr;
    
    % rotate and move base to xy = 0
    corbPlist = [];
    corbPlist2 = [];
    corbPlist3 = [];
    
    z_line = floor(mean(bPlist(:,3)));
    actinZ = [];
    
    for i = 1:size(filament_struct_expanded,2)
        
        actinZ = [actinZ;filament_struct_expanded(i).cor_filament_ext(:,3)];
        
    end
    
    actin_meanZ = mean(actinZ);

    
    
    for i = 1:size(bPlist,1)
    
        corbPlist(i,:) = actin_RotatePoints(bPlist(i,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),960,928,750);
        corbPlist2(i,:) = RotatePoints(corbPlist(i,:),0,0,-tomo_rot_psi,960,928,750);
        
        
        
        if actin_meanZ < z_line
            corbPlist3(i,:) = actin_RotatePoints(corbPlist2(i,:),pi,0,0,960,928,750);
        else
            corbPlist3(i,:) = corbPlist2(i,:);
        end
        
  

    end
    
    
    tomo_z_shift = round(mean(corbPlist3(:,3)));
%         [n,V,p] = affine_fit(corbPlist3);
%         [xx,yy] = meshgrid(0:1:1024,0:1:1024); 
%         d = -p*n; %'# dot product for less typing
%         z = (-n(1)*xx - n(2)*yy - d)./n(3);
%         surf(xx,yy,z,'FaceAlpha',0.5,'EdgeColor','none');
%         view(0,0);

    
    % rotate filaments
    for i = 1:size(filament_struct_expanded,2)
        cor = filament_struct_expanded(i).cor_filament_ext;
        psi_angle = filament_struct_expanded(i).psi_filament_ext;
        rot_xyz = filament_struct_expanded(i).rot_xyz_filament_ext;
        cor_rot = [];
        cor_rot2 = [];
        cor_rot3 = [];
        psi_angle_rot = [];
        rot_xyz_rot = [];
        rot_xyz_rot2 = [];
        rot_xyz_rot3 = [];
        for j = 1:size(cor,1)
            cor_rot(j,:) = actin_RotatePoints(cor(j,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),960,928,750); %in rad
            cor_rot2(j,:) = RotatePoints(cor_rot(j,:),0,0,-tomo_rot_psi,960,928,750); % in degree
            if actin_meanZ < z_line
                cor_rot3(j,:) = actin_RotatePoints(cor_rot2(j,:),pi,0,0,960,928,750);
            else
                cor_rot3(j,:) = cor_rot2(j,:);
            end
            
            if filament_struct_expanded(i).filament_processed_vec(j) ==1
                rot_xyz_rot(j,:) = actin_RotatePoints(rot_xyz(j,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),0,0,0);
                rot_xyz_rot2(j,:) = RotatePoints(rot_xyz_rot(j,:),0,0,-tomo_rot_psi,0,0,0);
                psi_angle_rot(j,1) = psi_angle(j,1)-rad2deg(tomo_rot_b(1,3))-tomo_rot_psi;
                if actin_meanZ < z_line
                    psi_angle_rot(j,1) = -(psi_angle_rot(j,1)-180);
                    rot_xyz_rot3(j,:) = actin_RotatePoints(rot_xyz_rot2(j,:),pi,0,0,0,0,0);
                else
                    rot_xyz_rot3(j,:) = rot_xyz_rot2(j,:);
                end
                
            else
                psi_angle_rot(j,1) = 0;
                rot_xyz_rot3(j,:) = rot_xyz(j,:);
            end
            
        end
        
        filament_struct_expanded(i).cor_filament_ext_psi_rot = cor_rot3;
        filament_struct_expanded(i).psi_filament_ext_rot = wrapTo180(psi_angle_rot);
        filament_struct_expanded(i).rot_xyz_filament_ext_rot = -rot_xyz_rot3 ;
    end
    
    %rotate stressF_vector
    
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
    
        for i = 1:size(filament_struct_expanded,2)
            
            cor = filament_struct_expanded(i).tomogram_stressF_vector;
            cor_rot = [];
            cor_rot2 = [];
            cor_rot3 = [];
            for j = 1:size(cor,1)
                cor_rot(j,:) = actin_RotatePoints(cor(j,:),tomo_rot_b(1,1),tomo_rot_b(1,2),tomo_rot_b(1,3),0,0,0); %in rad
                cor_rot2(j,:) = RotatePoints(cor_rot(j,:),0,0,-tomo_rot_psi,0,0,0); % in degree
                if actin_meanZ < z_line
                    cor_rot3(j,:) = actin_RotatePoints(cor_rot2(j,:),pi,0,0,0,0,0);
                else
                    cor_rot3(j,:) = cor_rot2(j,:);
                end
                
                
                
            end
            cor_rot3(:,3)=0; 
            filament_struct_expanded(i).tomogram_stressF_vector_rot = cor_rot3;
            
            
        end
    end
    
    %rotate stressF_vector to the side
    
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        if size(filament_struct_expanded(1).tomogram_stressF_vector,1) == 1
            n = filament_struct_expanded(1).tomogram_stressF_vector_rot;
            
        elseif size(filament_struct_expanded(1).tomogram_stressF_vector,1) > 1
            n = mean(filament_struct_expanded(1).tomogram_stressF_vector_rot);
        end
        
        
        n_norm = (n/norm(n));
        
        Zv = [1 0 0];
        corR = vrrotvec(n_norm,Zv);
        rotm = axang2rotm(corR);
        eul = rotm2eul(rotm,'XYZ');
        corV_angle = eul;
        
        
        %roate filament and again stressF_vector again
        
        for i = 1:size(filament_struct_expanded,2)
            cor = filament_struct_expanded(i).cor_filament_ext_psi_rot;
            psi_angle = filament_struct_expanded(i).psi_filament_ext_rot;
            rot_xyz = filament_struct_expanded(i).rot_xyz_filament_ext_rot;
            cor_rot = [];
            cor_rot2 = [];
            cor_rot3 = [];
            psi_angle_rot = [];
            rot_xyz_rot = [];
            rot_xyz_rot2 = [];
            rot_xyz_rot3 = [];
            for j = 1:size(cor,1)
                cor_rot(j,:) = actin_RotatePoints(cor(j,:),corV_angle(1,1),corV_angle(1,2),corV_angle(1,3),960,928,750); %in rad
                cor_rot2(j,:) = RotatePoints(cor_rot(j,:),0,0,-tomo_rot_psi,960,928,750); % in degree
                cor_rot3(j,:) = cor_rot2(j,:);
                
                if filament_struct_expanded(i).filament_processed_vec(j) ==1
                    rot_xyz_rot(j,:) = actin_RotatePoints(rot_xyz(j,:),corV_angle(1,1),corV_angle(1,2),corV_angle(1,3),0,0,0);
                    rot_xyz_rot2(j,:) = RotatePoints(rot_xyz_rot(j,:),0,0,-tomo_rot_psi,0,0,0);
                    psi_angle_rot(j,1) = psi_angle(j,1)-rad2deg(corV_angle(1,3))-tomo_rot_psi;
                    if ~ismember(k,[k])
                        psi_angle_rot(j,1) = -(psi_angle_rot(j,1)-180);
                        rot_xyz_rot3(j,:) = actin_RotatePoints(rot_xyz_rot2(j,:),pi,0,0,0,0,0);
                    else
                        rot_xyz_rot3(j,:) = rot_xyz_rot2(j,:);
                    end
                    
                else
                    psi_angle_rot(j,1) = 0;
                    rot_xyz_rot3(j,:) = rot_xyz(j,:);
                end
                
            end
            cor_rot3(:,3) = cor_rot3(:,3)-tomo_z_shift;
            filament_struct_expanded(i).cor_filament_ext_psi_rot_side = cor_rot3;
            filament_struct_expanded(i).psi_filament_ext_rot_side = wrapTo180(psi_angle_rot);
            filament_struct_expanded(i).rot_xyz_filament_ext_rot_side = rot_xyz_rot3 ;
        end
        
        
        %rotate stressF_vector
        
            
        for i = 1:size(filament_struct_expanded,2)
            
            cor = filament_struct_expanded(i).tomogram_stressF_vector_rot;
            cor_rot = [];
            cor_rot2 = [];
            cor_rot3 = [];
            for j = 1:size(cor,1)
                cor_rot(j,:) = actin_RotatePoints(cor(j,:),corV_angle(1,1),corV_angle(1,2),corV_angle(1,3),0,0,0); %in rad
                cor_rot2(j,:) = RotatePoints(cor_rot(j,:),0,0,-tomo_rot_psi,0,0,0); % in degree
                if ~ismember(k,[k])
                    cor_rot3(j,:) = actin_RotatePoints(cor_rot2(j,:),pi,0,0,0,0,0);
                else
                    cor_rot3(j,:) = cor_rot2(j,:);
                end
                
                
                
            end
            cor_rot3 = (cor_rot3/norm(cor_rot3));
            filament_struct_expanded(i).tomogram_stressF_vector_rot_side = cor_rot3;
            
            
        end
        
        
    end
    
    
    
    
 
    
%     % find the minZ in membrane and max Z of tomogram 
%     up = [];
%     down = [];
%     for i = 1:size(filament_struct_expanded,2)
%         up(i,:) = max(filament_struct_expanded(i).cor_filament_ext_psi_rot(:,3));
%     end
%     down = min(mPlist3(:,3));
%     upndown2 = [max(up) down max(up)-down];
%     for i = 1:size(filament_struct_expanded,2)
%         zzz = filament_struct_expanded(i).cor_filament_ext_psi_rot(:,3);
%         zzzR = (zzz-upndown2(2))./upndown2(3);
%         filament_struct_expanded(i).Z_thick = upndown2(3);
%         filament_struct_expanded(i).Z_Ration = zzzR;
%     end
    
    save(FilamentStruct2{k},'filament_struct_expanded');
    
    disp(k)

    

    
    
    
        
        
end
clear
    





%% plot check
actin_adhesion_base_path_cells_lessbent;

% load filaments
for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step5/filaments_160_tomo_' num2str(t) '.mat'];
end

plotfolder = '/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step5_plot/';
mkdir(plotfolder)
cd(plotfolder)
for k = [1:28 30:110 113:189]
    
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
    
        figure(k);
        hold on
        xlim([0 1024]);
        ylim([0 1024]);
        zlim([-100 512]);
        bbb = [];
        rrr = [];
        kkk = [];
        
        % plot filaments
        for i = 1:size(filament_struct_expanded,2)
            if filament_struct_expanded(i).filament_conf_com_score >= 0.6
                if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                    for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
                        if filament_struct_expanded(i).filament_spin_cat(j) == filament_struct_expanded(i).filament_spin_indx
                            bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(j,:)];
                            
                        elseif filament_struct_expanded(i).filament_spin_cat(j) == -(filament_struct_expanded(i).filament_spin_indx)
                            rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(j,:)];
                        elseif filament_struct_expanded(i).filament_spin_cat(j) == 0
                            kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(j,:)];
                        end
                    end
                    %               text(filament_struct_expanded(i).cor_filament_ext_psi_rot_side(j,1)...
                    %                ,filament_struct_expanded(i).cor_filament_ext_psi_rot_side(j,2)...
                    %                ,filament_struct_expanded(i).cor_filament_ext_psi_rot_side(j,3),num2str(i),'Color','green','FontSize',8)
                    
                    plot3(filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1)...
                        ,filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,2)...
                        ,filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,3),'k-')
                    
                    
                end
            end
        end
        
        if ~isempty(bbb)
            scatter3(bbb(:,1),bbb(:,2),bbb(:,3),40,'.','b')
        end
        if ~isempty(rrr)
            scatter3(rrr(:,1),rrr(:,2),rrr(:,3),40,'.','r')
        end
        if ~isempty(kkk)
            scatter3(kkk(:,1),kkk(:,2),kkk(:,3),40,'*','k')
        end
        
        
        
        

        
        set(gcf,'Position',[0 0 1024 1024])
        view(0,90)
        figname = ['tomo_' num2str(k)];
        saveas(gcf, fullfile(plotfolder,figname ),'png' )
        view(0,0)
        figname = ['tomoSide_' num2str(k)];
        saveas(gcf, fullfile(plotfolder,figname ),'png' )
        
        close(gcf)
    
    end

end

clear




