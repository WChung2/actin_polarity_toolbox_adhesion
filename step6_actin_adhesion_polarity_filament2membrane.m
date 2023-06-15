%% categorize actin to stress fiber orientation

% Load filaments and BasePath
actin_adhesion_base_path_cells_lessbent;

for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step5/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:189
    FilamentStruct2{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step6/filaments_160_tomo_' num2str(t) '.mat'];
end
mkdir /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step6/
cd /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step6/

pixelR = 0.220651*4;% conver to nm



for k = [1:28 30:110 113:189]
    % load filament structure
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        
        
        normalvector = filament_struct_expanded(1).tomogram_stressF_vector_rot_side(:,1:2);
        
        %calculate V_average
        for i = 1:size(filament_struct_expanded,2)
            avg_vector = [];
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                for j = [find(filament_struct_expanded(i).filament_processed_vec==1)']
                    
                    pointrot = [];
                    if filament_struct_expanded(i).filament_spin_cat(j)== filament_struct_expanded(i).filament_spin_indx
                        pointrot = filament_struct_expanded(i).rot_xyz_filament_ext_rot_side(j,:);
                    elseif filament_struct_expanded(i).filament_spin_cat(j)== -filament_struct_expanded(i).filament_spin_indx
                        pointrot = -filament_struct_expanded(i).rot_xyz_filament_ext_rot_side(j,:);
                    elseif (filament_struct_expanded(i).filament_spin_cat(j))/3 == filament_struct_expanded(i).filament_spin_indx
                        pointrot = filament_struct_expanded(i).rot_xyz_filament_ext_rot_side(j,:);
                    elseif (filament_struct_expanded(i).filament_spin_cat(j))/3 == -filament_struct_expanded(i).filament_spin_indx
                        pointrot = -filament_struct_expanded(i).rot_xyz_filament_ext_rot_side(j,:);
                    end
                   
                    
                    
                    avg_vector = [avg_vector;pointrot];
                end
                filament_struct_expanded(i).rot_xyz_filament_ext_rot_side_avg = mean(avg_vector);
            end
        end
        
        
        % calculate V_stress angle
        for i = 1:size(filament_struct_expanded,2)
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                for j = [find(filament_struct_expanded(i).filament_processed_vec==1)']
                    
                    pointrot = [];
                    if filament_struct_expanded(i).filament_spin_cat(j)== filament_struct_expanded(i).filament_spin_indx
                        pointrot = filament_struct_expanded(i).rot_xyz_filament_ext_rot_side(j,1:2);
                    elseif filament_struct_expanded(i).filament_spin_cat(j)== -filament_struct_expanded(i).filament_spin_indx
                        pointrot = -filament_struct_expanded(i).rot_xyz_filament_ext_rot_side(j,1:2);
                    else
                        pointrot = filament_struct_expanded(i).rot_xyz_filament_ext_rot_side(j,1:2);
                    end
                    delta_angle = [];
                    delta_angle_abs = [];
                    for q = 1: size(normalvector,1)
                        delta_angle(q) = acosd(dot(pointrot,normalvector(q,:))/(norm(pointrot)*norm(normalvector(q,:))));
                        delta_angle_abs(q) = acosd(abs(dot(pointrot,normalvector(q,:))/(norm(pointrot)*norm(normalvector(q,:)))));
                    end
                    filament_struct_expanded(i).point2V_stress_angle(j) = delta_angle((delta_angle_abs==min(delta_angle_abs)));
                    
                    
                end
            end
        end
        
        % categorze filaments 1==goining into -1==going away  to the distal end
        
        for i = 1:size(filament_struct_expanded,2)
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                filament_real_index = find(filament_struct_expanded(i).filament_processed_vec==1)';
                filamentrealindex = find(filament_struct_expanded(i).filament_processed_vec==1)';
                Fvector = filament_struct_expanded(i).cor_filament_ext_psi_rot_side(filament_real_index(end),:)-filament_struct_expanded(i).cor_filament_ext_psi_rot_side(filament_real_index(1),:);
                filament_struct_expanded(i).filament_avg_point2V_stress = rad2deg(atan2(mean(sind((filament_struct_expanded(i).point2V_stress_angle(filamentrealindex)))),...
                    mean(cosd((filament_struct_expanded(i).point2V_stress_angle(filamentrealindex))))));
                if abs(filament_struct_expanded(i).filament_avg_point2V_stress) > 45 && abs(filament_struct_expanded(i).filament_avg_point2V_stress) <135
                    filament_struct_expanded(i).filament_avg_move_cat= 2;
                elseif abs(filament_struct_expanded(i).filament_avg_point2V_stress) <= 45
                    filament_struct_expanded(i).filament_avg_move_cat= 1;
                elseif abs(filament_struct_expanded(i).filament_avg_point2V_stress) >= 135
                    filament_struct_expanded(i).filament_avg_move_cat= -1;
                end
                
                
                % define head and tail
                if filament_struct_expanded(i).filament_spin_cat(filamentrealindex(1)) == filament_struct_expanded(i).filament_spin_indx
                    pointrot  = filament_struct_expanded(i).rot_xyz_filament_ext_rot_side(filamentrealindex(1),:);
                    crossdeg = acosd(dot(pointrot,Fvector)/(norm(pointrot)*norm(Fvector)));
                    if crossdeg > 90
                        filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(1);
                        filament_struct_expanded(i).filamentPointindex = filamentrealindex(end);
                        
                    else
                        filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(end);
                        filament_struct_expanded(i).filamentPointindex = filamentrealindex(1);
                        
                    end
                elseif filament_struct_expanded(i).filament_spin_cat(filamentrealindex(1)) == -filament_struct_expanded(i).filament_spin_indx
                    pointrot  = -filament_struct_expanded(i).rot_xyz_filament_ext_rot(filamentrealindex(1),:);
                    crossdeg = acosd(dot(pointrot,Fvector)/(norm(pointrot)*norm(Fvector)));
                    if crossdeg > 90
                        filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(1);
                        filament_struct_expanded(i).filamentPointindex = filamentrealindex(end);
                    else
                        filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(end);
                        filament_struct_expanded(i).filamentPointindex = filamentrealindex(1);
                    end
                elseif filament_struct_expanded(i).filament_spin_cat(filamentrealindex(end)) == filament_struct_expanded(i).filament_spin_indx
                    pointrot  = filament_struct_expanded(i).rot_xyz_filament_ext_rot(filamentrealindex(end),:);
                    crossdeg = acosd(dot(pointrot,Fvector)/(norm(pointrot)*norm(Fvector)));
                    if crossdeg > 90
                        filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(1);
                        filament_struct_expanded(i).filamentPointindex = filamentrealindex(end);
                        
                    else
                        filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(end);
                        filament_struct_expanded(i).filamentPointindex = filamentrealindex(1);
                    end
                elseif filament_struct_expanded(i).filament_spin_cat(filamentrealindex(end)) == -filament_struct_expanded(i).filament_spin_indx
                    pointrot  = -filament_struct_expanded(i).rot_xyz_filament_ext_rot(filamentrealindex(1),:);
                    crossdeg = acosd(dot(pointrot,Fvector)/(norm(pointrot)*norm(Fvector)));
                    if crossdeg > 90
                        filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(1);
                        filament_struct_expanded(i).filamentPointindex = filamentrealindex(end);
                    else
                        filament_struct_expanded(i).filamentBarbedindex = filamentrealindex(end);
                        filament_struct_expanded(i).filamentPointindex = filamentrealindex(1);
                    end
                    
                end
                
            end
        end
    end
    
    
    
    save(FilamentStruct2{k},'filament_struct_expanded');
    disp(k)
end
clear



%% plot avg angle
actin_adhesion_base_path_cells_lessbent;

% load filaments
for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step6/filaments_160_tomo_' num2str(t) '.mat'];
end

plotfolder = '/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step6_plot/';

mkdir(plotfolder)
cd(plotfolder)

filament_com_conf_thres = 0.6;
cmap = turbo(360);

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
        ggg = [];
        kkk = [];
        qqq = [];
        
        % plot filaments in 3 colors
        for i = 1:size(filament_struct_expanded,2)
            
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
                    if filament_struct_expanded(i).filament_avg_move_cat == 1
                        bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(filament_struct_expanded(i).filament_processed_vec==1,:)];
                        
                    elseif filament_struct_expanded(i).filament_avg_move_cat == -1
                        rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(filament_struct_expanded(i).filament_processed_vec==1,:)];
                        
                    elseif filament_struct_expanded(i).filament_avg_move_cat == 2
                        ggg = [ggg; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(filament_struct_expanded(i).filament_processed_vec==1,:)];
                        
                        
                    end
                    kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(logical(filament_struct_expanded(i).tropoIndex),:)];
                    qqq = [qqq; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(filament_struct_expanded(i).filament_processed_vec==0,:)];
                    
                    
                    
%                     text(filament_struct_expanded(i).cor_filament_ext_psi_rot_side(end,1)...
%                         ,filament_struct_expanded(i).cor_filament_ext_psi_rot_side(end,2)...
%                         ,filament_struct_expanded(i).cor_filament_ext_psi_rot_side(end,3),num2str(i))
                    
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
        if ~isempty(ggg)
            scatter3(ggg(:,1),ggg(:,2),ggg(:,3),40,'.','g')
        end
        if ~isempty(kkk)
            scatter3(kkk(:,1),kkk(:,2),kkk(:,3),70,'o','m','filled')
        end
%         if ~isempty(qqq)
%             scatter3(qqq(:,1),qqq(:,2),qqq(:,3),20,'.','k')
%         end        
        % plot filaments in 360 colors
        %     for i = 1:size(filament_struct_expanded,2)
        %
        %         if ~isempty(filament_struct_expanded(i).filament_spin_indx)
        % %             if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
        %                 for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
        %                     color = round(wrapTo360(filament_struct_expanded(i).point2membrane_angle(j)));
        %                   if filament_struct_expanded(i).filament_move_cat(j) ~= 0
        %                       scatter3(filament_struct_expanded(i).cor_filament_ext_psi_rot(j,1)...
        %                           ,filament_struct_expanded(i).cor_filament_ext_psi_rot(j,2)...
        %                           ,filament_struct_expanded(i).cor_filament_ext_psi_rot(j,3),50,cmap(color,:),'.')
        %
        %
        %                   elseif filament_struct_expanded(i).filament_move_cat(j) == 0
        %                       scatter3(filament_struct_expanded(i).cor_filament_ext_psi_rot(j,1)...
        %                           ,filament_struct_expanded(i).cor_filament_ext_psi_rot(j,2)...
        %                           ,filament_struct_expanded(i).cor_filament_ext_psi_rot(j,3),'*','k')
        %                   end
        %
        %                 end
        %
        %                   plot3(filament_struct_expanded(i).cor_filament_ext_psi_rot(:,1)...
        %                       ,filament_struct_expanded(i).cor_filament_ext_psi_rot(:,2)...
        %                       ,filament_struct_expanded(i).cor_filament_ext_psi_rot(:,3),'k-')
        %
        % %             end
        %
        %
        %
        %          end
        %     end
        
        
        set(gcf,'Position',[0 0 1024 1024])
        view(0,90)
        figname = ['tomo_new' num2str(k)];
        saveas(gcf, fullfile(plotfolder,figname ),'png' )
        view(0,0)
        figname = ['tomoSide_new' num2str(k)];
        saveas(gcf, fullfile(plotfolder,figname ),'png' )
        
        close(gcf)
    end

end

clear


%% subplot to check in 2D

actin_adhesion_base_path_cells_lessbent;

% load filaments
for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step6/filaments_160_tomo_' num2str(t) '.mat'];
end

plotfolder = '/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step6_plot';

mkdir(plotfolder)
cd(plotfolder)

filament_com_conf_thres = 0.6;
cmap = turbo(360);


wt_good = [6:31 33 34 96:144 146:150];
v_good = [74:95];
zyx_good = [35:38 40:45 151:189];
double = [1 2 4 5 46:51 52 54:64 66:73];


figure(1);
set(gcf,'Position',[1 1 4096 1024])
counter = 1;
% top view
for k = wt_good
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        subplot(7,12,counter)
        hold on
        xlim([0 1024]);
        ylim([0 1024]);
        bbb = [];
        rrr = [];
        ggg = [];
        kkk = [];

        
        
        % plot filaments in 3 colors
        
        for i = 1:size(filament_struct_expanded,2)
            
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
                    for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
                        if filament_struct_expanded(i).filament_avg_move_cat == 1
                            bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == -1
                            rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == 2
                            ggg = [ggg; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        end
                        kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(logical(filament_struct_expanded(i).tropoIndex),1:2)];
                    end
                    
                    
                end
                
                
                
            end
        end
        if ~isempty(bbb)
            scatter(bbb(:,1),bbb(:,2),7,'.','b')
        end
        if ~isempty(rrr)
            scatter(rrr(:,1),rrr(:,2),7,'.','r')
        end
        if ~isempty(ggg)
            scatter(ggg(:,1),ggg(:,2),7,'.','g')
        end
        if ~isempty(kkk)
            scatter(kkk(:,1),kkk(:,2),3,'o','m')
        end
        
        
        
        title(sprintf('tomo_ %d',k));
        disp(k)
        counter = counter+1;
    end

end
figname = ('all_wt_tropo');
saveas(gcf, fullfile(plotfolder,figname ),'png' )
close(gcf)

figure(1);
set(gcf,'Position',[1 1 4096 1024])
counter = 1;
% side view
for k = wt_good
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        subplot(7,12,counter)
        hold on
        xlim([0 1024]);
        ylim([-10 400]);
        bbb = [];
        rrr = [];
        ggg = [];
        kkk = [];

        
        
        % plot filaments in 3 colors
        
        for i = 1:size(filament_struct_expanded,2)
            
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
                    for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
                        if filament_struct_expanded(i).filament_avg_move_cat == 1
                            bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == -1
                            rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == 2
                            ggg = [ggg; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        end
                        kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(logical(filament_struct_expanded(i).tropoIndex),[1 3])];
                        
                    end
                    
                    
                end
                
                
                
            end
        end
        if ~isempty(bbb)
            scatter(bbb(:,1),bbb(:,2),7,'.','b')
        end
        if ~isempty(rrr)
            scatter(rrr(:,1),rrr(:,2),7,'.','r')
        end
        if ~isempty(ggg)
            scatter(ggg(:,1),ggg(:,2),7,'.','g')
        end
        if ~isempty(kkk)
            scatter(kkk(:,1),kkk(:,2),3,'o','m')
        end
        
        
        
        title(sprintf('tomo_ %d',k));
        disp(k)
        counter = counter+1;
    end

end
figname = ('all_wt_tropo_side');
saveas(gcf, fullfile(plotfolder,figname ),'png' )
close(gcf)


figure(1);
set(gcf,'Position',[1 1 4096 1024])
counter = 1;
% top view
for k = v_good
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        subplot(7,12,counter)
        hold on
        xlim([0 1024]);
        ylim([0 1024]);
        bbb = [];
        rrr = [];
        ggg = [];
        kkk = [];

        
        
        % plot filaments in 3 colors
        
        for i = 1:size(filament_struct_expanded,2)
            
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
                    for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
                        if filament_struct_expanded(i).filament_avg_move_cat == 1
                            bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == -1
                            rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == 2
                            ggg = [ggg; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        end
                        kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(logical(filament_struct_expanded(i).tropoIndex),1:2)];
                    end
                    
                    
                end
                
                
                
            end
        end
        if ~isempty(bbb)
            scatter(bbb(:,1),bbb(:,2),7,'.','b')
        end
        if ~isempty(rrr)
            scatter(rrr(:,1),rrr(:,2),7,'.','r')
        end
        if ~isempty(ggg)
            scatter(ggg(:,1),ggg(:,2),7,'.','g')
        end
        if ~isempty(kkk)
            scatter(kkk(:,1),kkk(:,2),3,'o','m')
        end
        
        
        
        title(sprintf('tomo_ %d',k));
        disp(k)
        counter = counter+1;
    end

end
figname = ('all_vasp_tropo');
saveas(gcf, fullfile(plotfolder,figname ),'png' )
close(gcf)

figure(1);
set(gcf,'Position',[1 1 4096 1024])
counter = 1;
% side view
for k = v_good
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        subplot(7,12,counter)
        hold on
        xlim([0 1024]);
        ylim([-10 400]);
        bbb = [];
        rrr = [];
        ggg = [];
        kkk = [];

        
        
        % plot filaments in 3 colors
        
        for i = 1:size(filament_struct_expanded,2)
            
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
                    for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
                        if filament_struct_expanded(i).filament_avg_move_cat == 1
                            bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == -1
                            rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == 2
                            ggg = [ggg; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        end
                        kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(logical(filament_struct_expanded(i).tropoIndex),[1 3])];
                        
                    end
                    
                    
                end
                
                
                
            end
        end
        if ~isempty(bbb)
            scatter(bbb(:,1),bbb(:,2),7,'.','b')
        end
        if ~isempty(rrr)
            scatter(rrr(:,1),rrr(:,2),7,'.','r')
        end
        if ~isempty(ggg)
            scatter(ggg(:,1),ggg(:,2),7,'.','g')
        end
        if ~isempty(kkk)
            scatter(kkk(:,1),kkk(:,2),3,'o','m')
        end
        
        
        
        title(sprintf('tomo_ %d',k));
        disp(k)
        counter = counter+1;
    end

end
figname = ('all_vasp_tropo_side');
saveas(gcf, fullfile(plotfolder,figname ),'png' )
close(gcf)

figure(1);
set(gcf,'Position',[1 1 4096 1024])
counter = 1;
% top view
for k = zyx_good
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        subplot(7,12,counter)
        hold on
        xlim([0 1024]);
        ylim([0 1024]);
        bbb = [];
        rrr = [];
        ggg = [];
        kkk = [];

        
        
        % plot filaments in 3 colors
        
        for i = 1:size(filament_struct_expanded,2)
            
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
                    for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
                        if filament_struct_expanded(i).filament_avg_move_cat == 1
                            bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == -1
                            rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == 2
                            ggg = [ggg; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        end
                        kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(logical(filament_struct_expanded(i).tropoIndex),1:2)];
                    end
                    
                    
                end
                
                
                
            end
        end
        if ~isempty(bbb)
            scatter(bbb(:,1),bbb(:,2),7,'.','b')
        end
        if ~isempty(rrr)
            scatter(rrr(:,1),rrr(:,2),7,'.','r')
        end
        if ~isempty(ggg)
            scatter(ggg(:,1),ggg(:,2),7,'.','g')
        end
        if ~isempty(kkk)
            scatter(kkk(:,1),kkk(:,2),3,'o','m')
        end
        
        
        
        title(sprintf('tomo_ %d',k));
        disp(k)
        counter = counter+1;
    end

end
figname = ('all_zyx_tropo');
saveas(gcf, fullfile(plotfolder,figname ),'png' )
close(gcf)

figure(1);
set(gcf,'Position',[1 1 4096 1024])
counter = 1;
% side view
for k = zyx_good
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        subplot(7,12,counter)
        hold on
        xlim([0 1024]);
        ylim([-10 400]);
        bbb = [];
        rrr = [];
        ggg = [];
        kkk = [];

        
        
        % plot filaments in 3 colors
        
        for i = 1:size(filament_struct_expanded,2)
            
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
                    for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
                        if filament_struct_expanded(i).filament_avg_move_cat == 1
                            bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == -1
                            rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == 2
                            ggg = [ggg; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        end
                        kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(logical(filament_struct_expanded(i).tropoIndex),[1 3])];
                        
                    end
                    
                    
                end
                
                
                
            end
        end
        if ~isempty(bbb)
            scatter(bbb(:,1),bbb(:,2),7,'.','b')
        end
        if ~isempty(rrr)
            scatter(rrr(:,1),rrr(:,2),7,'.','r')
        end
        if ~isempty(ggg)
            scatter(ggg(:,1),ggg(:,2),7,'.','g')
        end
        if ~isempty(kkk)
            scatter(kkk(:,1),kkk(:,2),3,'o','m')
        end
        
        
        
        title(sprintf('tomo_ %d',k));
        disp(k)
        counter = counter+1;
    end

end
figname = ('all_zyx_tropo_side');
saveas(gcf, fullfile(plotfolder,figname ),'png' )
close(gcf)

figure(1);
set(gcf,'Position',[1 1 4096 1024])
counter = 1;
% top view
for k = double
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        subplot(7,12,counter)
        hold on
        xlim([0 1024]);
        ylim([0 1024]);
        bbb = [];
        rrr = [];
        ggg = [];
        kkk = [];

        
        
        % plot filaments in 3 colors
        
        for i = 1:size(filament_struct_expanded,2)
            
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
                    for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
                        if filament_struct_expanded(i).filament_avg_move_cat == 1
                            bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == -1
                            rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == 2
                            ggg = [ggg; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,1:2)];
                            
                        end
                        kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(logical(filament_struct_expanded(i).tropoIndex),1:2)];
                    end
                    
                    
                end
                
                
                
            end
        end
        if ~isempty(bbb)
            scatter(bbb(:,1),bbb(:,2),7,'.','b')
        end
        if ~isempty(rrr)
            scatter(rrr(:,1),rrr(:,2),7,'.','r')
        end
        if ~isempty(ggg)
            scatter(ggg(:,1),ggg(:,2),7,'.','g')
        end
        if ~isempty(kkk)
            scatter(kkk(:,1),kkk(:,2),3,'o','m')
        end
        
        
        
        title(sprintf('tomo_ %d',k));
        disp(k)
        counter = counter+1;
    end

end
figname = ('all_double_tropo');
saveas(gcf, fullfile(plotfolder,figname ),'png' )
close(gcf)

figure(1);
set(gcf,'Position',[1 1 4096 1024])
counter = 1;
% side view
for k = double
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        subplot(7,12,counter)
        hold on
        xlim([0 1024]);
        ylim([-10 400]);
        bbb = [];
        rrr = [];
        ggg = [];
        kkk = [];

        
        
        % plot filaments in 3 colors
        
        for i = 1:size(filament_struct_expanded,2)
            
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
                    for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
                        if filament_struct_expanded(i).filament_avg_move_cat == 1
                            bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == -1
                            rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        elseif filament_struct_expanded(i).filament_avg_move_cat == 2
                            ggg = [ggg; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(:,[1 3])];
                            
                        end
                        kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side(logical(filament_struct_expanded(i).tropoIndex),[1 3])];
                        
                    end
                    
                    
                end
                
                
                
            end
        end
        if ~isempty(bbb)
            scatter(bbb(:,1),bbb(:,2),7,'.','b')
        end
        if ~isempty(rrr)
            scatter(rrr(:,1),rrr(:,2),7,'.','r')
        end
        if ~isempty(ggg)
            scatter(ggg(:,1),ggg(:,2),7,'.','g')
        end
        if ~isempty(kkk)
            scatter(kkk(:,1),kkk(:,2),3,'o','m')
        end
        
        
        
        title(sprintf('tomo_ %d',k));
        disp(k)
        counter = counter+1;
    end

end
figname = ('all_double_tropo_side');
saveas(gcf, fullfile(plotfolder,figname ),'png' )
close(gcf)



clear


%% debug plot

actin_adhesion_base_path_cells_lessbent;

% load filaments
for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step6/filaments_160_tomo_' num2str(t) '.mat'];
end

plotfolder = '/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step6_plot/';

mkdir(plotfolder)
cd(plotfolder)

filament_com_conf_thres = 0.6;

for k = 1:189
    load(FilamentStruct{k},'filament_struct_expanded');
    if ~isempty(filament_struct_expanded(1).tomogram_stressF_vector)
        
        figure(k);
        hold on
        xlim([0 1024]);
        ylim([0 1024]);
        zlim([-100 512]);
        bbb = [];
        rrr = [];
        ggg = [];
        kkk = [];
        vbb = [];
        vrr = [];
        % plot filaments in 3 colors
        for i = 1:size(filament_struct_expanded,2)
            
            if ~isempty(filament_struct_expanded(i).filament_spin_indx)
                if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
                    if filament_struct_expanded(i).filament_avg_move_cat == 1
                        bbb = [bbb; filament_struct_expanded(i).cor_filament_ext_psi_rot_side];
                        vbb = [vbb;(filament_struct_expanded(i).rot_xyz_filament_ext_rot_side)];
                    elseif filament_struct_expanded(i).filament_avg_move_cat == -1
                        rrr = [rrr; filament_struct_expanded(i).cor_filament_ext_psi_rot_side];
                        vrr = [vrr;(filament_struct_expanded(i).rot_xyz_filament_ext_rot_side)];
                    elseif filament_struct_expanded(i).filament_avg_move_cat == 2
                        ggg = [ggg; filament_struct_expanded(i).cor_filament_ext_psi_rot_side];
                        
                    elseif filament_struct_expanded(i).filament_avg_move_cat == 0
                        kkk = [kkk; filament_struct_expanded(i).cor_filament_ext_psi_rot_side];
                        
                    end
                    
                    
                    
                    
%                     text(filament_struct_expanded(i).cor_filament_ext_psi_rot_side(end,1)...
%                         ,filament_struct_expanded(i).cor_filament_ext_psi_rot_side(end,2)...
%                         ,filament_struct_expanded(i).cor_filament_ext_psi_rot_side(end,3),num2str(i))
                    
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
        if ~isempty(ggg)
            scatter3(ggg(:,1),ggg(:,2),ggg(:,3),40,'.','g')
        end
        if ~isempty(kkk)
            scatter3(kkk(:,1),kkk(:,2),kkk(:,3),40,'*','k')
        end
        quiver3(bbb(:,1),bbb(:,2),bbb(:,3),vbb(:,1),vbb(:,2),vbb(:,3))
        quiver3(rrr(:,1),rrr(:,2),rrr(:,3),vrr(:,1),vrr(:,2),vrr(:,3))
        % plot filaments in 360 colors
        %     for i = 1:size(filament_struct_expanded,2)
        %
        %         if ~isempty(filament_struct_expanded(i).filament_spin_indx)
        % %             if filament_struct_expanded(i).filament_conf_com_score > filament_com_conf_thres
        %                 for j = 1: size(filament_struct_expanded(i).filament_spin_cat,2)
        %                     color = round(wrapTo360(filament_struct_expanded(i).point2membrane_angle(j)));
        %                   if filament_struct_expanded(i).filament_move_cat(j) ~= 0
        %                       scatter3(filament_struct_expanded(i).cor_filament_ext_psi_rot(j,1)...
        %                           ,filament_struct_expanded(i).cor_filament_ext_psi_rot(j,2)...
        %                           ,filament_struct_expanded(i).cor_filament_ext_psi_rot(j,3),50,cmap(color,:),'.')
        %
        %
        %                   elseif filament_struct_expanded(i).filament_move_cat(j) == 0
        %                       scatter3(filament_struct_expanded(i).cor_filament_ext_psi_rot(j,1)...
        %                           ,filament_struct_expanded(i).cor_filament_ext_psi_rot(j,2)...
        %                           ,filament_struct_expanded(i).cor_filament_ext_psi_rot(j,3),'*','k')
        %                   end
        %
        %                 end
        %
        %                   plot3(filament_struct_expanded(i).cor_filament_ext_psi_rot(:,1)...
        %                       ,filament_struct_expanded(i).cor_filament_ext_psi_rot(:,2)...
        %                       ,filament_struct_expanded(i).cor_filament_ext_psi_rot(:,3),'k-')
        %
        % %             end
        %
        %
        %
        %          end
        %     end
        
        
        set(gcf,'Position',[0 0 1024 1024])
        view(0,90)
        figname = ['tomo_' num2str(k)];
%         saveas(gcf, fullfile(plotfolder,figname ),'png' )
        view(0,0)
        figname = ['tomoSide_' num2str(k)];
%         saveas(gcf, fullfile(plotfolder,figname ),'png' )
        
        close(gcf)
    end

end

clear
