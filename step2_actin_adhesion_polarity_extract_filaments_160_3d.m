%% Build filaments structure 


% Load actin data matrix
load('/home/Medalia/Projects7/WChung/actin_adhesion/mapping/actin_polarity_mapping_160_step_1_2.mat','actin_data_matrix');

%load tropomyocin index
tropoMatirx = importdata('/home/Medalia/Projects7/WChung/actin_adhesion/relion4_apt/Refine3D/job080/run_data_mod.star');
tropoIndex = zeros(size(tropoMatirx.textdata,1),1);
for k = 1:size(tropoMatirx.textdata,1)
    index = [];
    index = str2double(tropoMatirx.textdata{k,1}(1:strfind(tropoMatirx.textdata{k,1},'@')-1));
    tropoIndex(k,1) = index;

end

% Define filament length in terms of segments to analyze threshold
filament_length_thres = 3;

% Define filament processed fraction threshold
% filament_processed_frac_thres = 0.5;

for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1

% Process bundles
for k=1:189

    % Select only those segments that were used in the respective tomogram and in Relion 3D refine
    indx_segments = find(actin_data_matrix(:,1) == k); %& actin_data_matrix(:,4) == 1);
    
    % Reduce actin data matrix
    actin_data_matrix_red = actin_data_matrix(indx_segments,:);
    
    % input tropo
    for i = 1:size(actin_data_matrix_red,1)
        if ismember(actin_data_matrix_red(i,3),tropoIndex)
            actin_data_matrix_red(i,39) = 1;
        else
            actin_data_matrix_red(i,39) = 0;
        end
        
    end

    
    % Extract filament assignment of each segment in the respective tomogram
    filament_vect = actin_data_matrix_red(:,38);
    
    % Extract filaments and filament distributions
    filament_struct = [];
    rot_distribution = [];
    ccc_distribution = [];
    tomo_distribution = [];
    segment_sensitivity_distribution = [];
    F_tropoIndex = [];
    zaehler_f = 1;
        
    
    
    for i=1:max(filament_vect(:))
        
        % Extract filament indices
        indx_filament_ext = find(filament_vect==i);
        
        % Extract real_processed_filament indices
        
        indx_filament_ext_real = find(actin_data_matrix_red(indx_filament_ext,4) == 1);
                        
        % Measure filament length
        filament_length = size(indx_filament_ext,1);
        
        % Measure filament processed fraction
        %filament_processed_fraction = size(indx_filament_ext_real,1)/filament_length;
        
        if ~isempty(indx_filament_ext) && size(indx_filament_ext_real,1) >= filament_length_thres %&& (filament_processed_fraction > filament_processed_frac_thres)
        
              filament_struct(zaehler_f).filament_indx = zaehler_f;
              
              % Extract coordinates of alignment corrected segments
              cor_filament_ext = [actin_data_matrix_red(indx_filament_ext,8)./4 actin_data_matrix_red(indx_filament_ext,9)./4 actin_data_matrix_red(indx_filament_ext,10)./4];
              cor_filament_delta = [actin_data_matrix_red(indx_filament_ext,14)./4 actin_data_matrix_red(indx_filament_ext,15)./4 actin_data_matrix_red(indx_filament_ext,16)./4];
              cor_filament_inv_delta = [actin_data_matrix_red(indx_filament_ext,21)./4 actin_data_matrix_red(indx_filament_ext,22)./4 actin_data_matrix_red(indx_filament_ext,23)./4];
              % Extract actin directions
              phi_filament_ext = actin_data_matrix_red(indx_filament_ext,18);
              psi_filament_ext = actin_data_matrix_red(indx_filament_ext,19);
              theta_filament_ext = actin_data_matrix_red(indx_filament_ext,20);
              
              % Extract actin directions
              rot_filament_ext = actin_data_matrix_red(indx_filament_ext,24:26);
              rot_distribution = [rot_distribution; rot_filament_ext];
              
              % Extract ccc distribution
              ccc_filament_ext = actin_data_matrix_red(indx_filament_ext,35);
              ccc_distribution = [ccc_distribution; ccc_filament_ext];
              
              % Extract tomo distribution
              tomo_filament_ext = actin_data_matrix_red(indx_filament_ext,1);
              tomo_distribution = [tomo_distribution; tomo_filament_ext];
              
              % Extract filament sensitivity distribution
              segment_sensitivity_filament_ext = actin_data_matrix_red(indx_filament_ext,37);
              segment_sensitivity_distribution = [segment_sensitivity_distribution; segment_sensitivity_filament_ext];
              filament_processed_vec = actin_data_matrix_red(indx_filament_ext,4);
              
              % Extract tropoIndex
              F_tropoIndex = actin_data_matrix_red(indx_filament_ext,39);
              
              % Extract actin vector
              vect_filament_ext = [actin_data_matrix_red(indx_filament_ext,24) actin_data_matrix_red(indx_filament_ext,25) actin_data_matrix_red(indx_filament_ext,26)];
              
              filament_struct(zaehler_f).cor_filament_ext = cor_filament_ext;
              filament_struct(zaehler_f).cor_filament_delta = cor_filament_delta;
              filament_struct(zaehler_f).cor_filament_inv_delta = cor_filament_inv_delta;
              filament_struct(zaehler_f).phi_filament_ext = phi_filament_ext;
              filament_struct(zaehler_f).psi_filament_ext = psi_filament_ext;
              filament_struct(zaehler_f).theta_filament_ext = theta_filament_ext;
              filament_struct(zaehler_f).rot_xyz_filament_ext = rot_filament_ext;
              filament_struct(zaehler_f).ccc_filament_ext = ccc_filament_ext;
              filament_struct(zaehler_f).segment_sensitivity_filament_ext = segment_sensitivity_filament_ext;
              filament_struct(zaehler_f).vect_filament_ext = vect_filament_ext;
              filament_struct(zaehler_f).filament_processed_vec = filament_processed_vec;
              %filament_struct(zaehler_f).filament_processed_fraction = filament_processed_fraction;
              filament_struct(zaehler_f).filament_length = size(indx_filament_ext,1);
              filament_struct(zaehler_f).tropoIndex = F_tropoIndex;
              zaehler_f = zaehler_f + 1;
              
        end
        
    end
    
   
    
    save(FilamentStruct{k},'filament_struct');
    
    disp(k)

end
clear
    

    
    
%% test plot 

for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1_plot/
cd /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1_plot/


for k=1:189
    
    load(FilamentStruct{k},'filament_struct');
    
    figure(k)
    hold on
    xlim([0 1024])
    ylim([0 1024])
    zlim([0 1024])
    rrr = [];
    kkk = [];
    for i =1:size(filament_struct,2)
        filament_process = [];
        filament_process = filament_struct(i).filament_processed_vec;
        rrr = [rrr;filament_struct(i).cor_filament_ext(filament_process==1,:)];
        kkk = [kkk;filament_struct(i).cor_filament_ext(filament_process==0,:)];
        plot3(filament_struct(i).cor_filament_ext(:,1)...
                    ,filament_struct(i).cor_filament_ext(:,2)...
                    ,filament_struct(i).cor_filament_ext(:,3),'b')
        text(filament_struct(i).cor_filament_ext(1,1)...
           ,filament_struct(i).cor_filament_ext(1,2)...
           ,filament_struct(i).cor_filament_ext(1,3),num2str(i),'Color','green','FontSize',8)
       
       
        
    end
    scatter3(rrr(:,1),rrr(:,2),rrr(:,3),'.','r')
    scatter3(kkk(:,1),kkk(:,2),kkk(:,3),'*','k')

    figname = ['tomo_' num2str(k)];
    set(gcf,'Position',[1 1 1024 1024],'Renderer','painters')
    saveas(gca, fullfile('/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1_plot/',figname ),'png' )
    view(0,0)
    figname = ['tomo_side' num2str(k)];
    saveas(gca, fullfile('/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1_plot/',figname ),'png' )

    close(gcf)
    
end

clear


%% plot orientation


for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1_plot_psi_colored/
cd /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1_plot_psi_colored/


for k=1:189
    
    load(FilamentStruct{k},'filament_struct');
    
    figure(k)
    hold on
    xlim([0 1024])
    ylim([0 1024])
    zlim([0 1024])
    turbo_wrap = turbo(360);
    colormap(turbo_wrap)
    for i =1:size(filament_struct,2)
        filament_cor = [];
        filament_process = [];
        filament_cor = filament_struct(i).cor_filament_ext;
        filament_process = filament_struct(i).filament_processed_vec;
        psi = round(filament_struct(i).psi_filament_ext(filament_process==1)+180);
        [idx,C] = kmeans(psi,2,'Replicates',20);
        if size(find(idx==1),1) >= size(find(idx==2),1)
            psi(idx==2) = wrapTo360(psi(idx==2)-180);
        else
            psi(idx==1) = wrapTo360(psi(idx==1)-180);
        end
        psi(psi==0)=360;
        scatter3(filament_struct(i).cor_filament_ext(filament_process==1,1),...
            filament_struct(i).cor_filament_ext(filament_process==1,2),...
            filament_struct(i).cor_filament_ext(filament_process==1,3),15,turbo_wrap(psi,:),'filled')
        
%         scatter3(filament_struct(i).cor_filament_ext(filament_process==0,1),...
%             filament_struct(i).cor_filament_ext(filament_process==0,2),...
%             filament_struct(i).cor_filament_ext(filament_process==0,3),'*','k')
        plot3(filament_struct(i).cor_filament_ext(:,1)...
            ,filament_struct(i).cor_filament_ext(:,2)...
            ,filament_struct(i).cor_filament_ext(:,3),'k')
%         text(filament_struct(i).cor_filament_ext(1,1)...
%             ,filament_struct(i).cor_filament_ext(1,2)...
%             ,filament_struct(i).cor_filament_ext(1,3),num2str(i),'Color','green','FontSize',8)
%         
        
        
    end
    colorbar;
    colorbar('Ticks',[0:0.25:1],'TickLabels',{'0','90','180','270','360'})
    figname = ['tomo_' num2str(k)];
    set(gcf,'Position',[1 1 1024 1024],'Renderer','painters')
    saveas(gca, fullfile('/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1_plot_psi_colored/',figname ),'png' )
    close(gcf)
    
end

clear

        
%% find the longest good filaments

for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step1/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:189
    FilamentStruct2{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step2/filaments_160_tomo_' num2str(t) '.mat'];
end

boxsize = 40*(0.220651*4);%in nm
Distance_next = 200000*boxsize;
filament_length_thres = 3;


for k=1:189
    zaehler_f = 1;
    load(FilamentStruct{k},'filament_struct');
    filament_struct_expanded = [];
    for i= 1:size(filament_struct,2)
        cor = [];
        realcor = [];
        cor = filament_struct(i).cor_filament_ext;
        filament_realindx = find(filament_struct(i).filament_processed_vec==1);
        realcor = cor(filament_realindx,:);
        pDist = [];
        for j = 1:size(realcor,1)-1
            pDist(j,1) = pdist2(realcor(j,:),realcor(j+1,:));
            if pDist(j,1) > Distance_next
                pDist(j,2) = -1;
            else
                pDist(j,2) = 0;
                
            end
        end
        breakpoint = filament_realindx(pDist(:,2)==(-1));
        startpoint = 0;
        ultrafilament = {};
        counter = 1;
        if ~isempty(breakpoint)
            for j = 1:size(breakpoint,1)
                ultrafilament{counter}= (startpoint+1):1:(breakpoint(j));
                startpoint = breakpoint(j);
                counter = counter+1;
            end
        end
        ultrafilament{counter} = (startpoint+1):1:size(filament_struct(i).cor_filament_ext,1);
        
        for j =1: size(ultrafilament,2)
%             filament_processed_fraction =...
%                 size(find(filament_struct(i).filament_processed_vec(ultrafilament{j},:)==1),1)/...
%                 size(filament_struct(i).filament_processed_vec(ultrafilament{j},:),1);
            if size(find(filament_struct(i).filament_processed_vec(ultrafilament{j},:)),1)>= filament_length_thres  %&& (filament_processed_fraction >= filament_processed_frac_thres)
                
                filament_struct_expanded(zaehler_f).filament_indx = zaehler_f;
                filament_struct_expanded(zaehler_f).cor_filament_ext = filament_struct(i).cor_filament_ext(ultrafilament{j},:);
                filament_struct_expanded(zaehler_f).cor_filament_delta = filament_struct(i).cor_filament_delta(ultrafilament{j},:);
                filament_struct_expanded(zaehler_f).cor_filament_inv_delta = filament_struct(i).cor_filament_inv_delta(ultrafilament{j},:);
                filament_struct_expanded(zaehler_f).phi_filament_ext = filament_struct(i).phi_filament_ext(ultrafilament{j},:);
                filament_struct_expanded(zaehler_f).psi_filament_ext = filament_struct(i).psi_filament_ext(ultrafilament{j},:);
                filament_struct_expanded(zaehler_f).theta_filament_ext = filament_struct(i).theta_filament_ext(ultrafilament{j},:);
                filament_struct_expanded(zaehler_f).rot_xyz_filament_ext = filament_struct(i).rot_xyz_filament_ext(ultrafilament{j},:);
                filament_struct_expanded(zaehler_f).ccc_filament_ext = filament_struct(i).ccc_filament_ext(ultrafilament{j},:);
                filament_struct_expanded(zaehler_f).segment_sensitivity_filament_ext = filament_struct(i).segment_sensitivity_filament_ext(ultrafilament{j},:);
                filament_struct_expanded(zaehler_f).vect_filament_ext = filament_struct(i).vect_filament_ext(ultrafilament{j},:);
                filament_struct_expanded(zaehler_f).filament_processed_vec = filament_struct(i).filament_processed_vec(ultrafilament{j},:);
                %filament_struct_expanded(zaehler_f).filament_processed_fraction = filament_processed_fraction;
                filament_struct_expanded(zaehler_f).filament_length = size(filament_struct(i).filament_processed_vec(ultrafilament{j},:),1);
                filament_struct_expanded(zaehler_f).tropoIndex = filament_struct(i).tropoIndex(ultrafilament{j},:);
                zaehler_f = zaehler_f + 1;
            end
            
        end
        
        
    end
    
    
 
    save(FilamentStruct2{k},'filament_struct_expanded');
    disp(k)
    
end
clear
%% find regression line for each filament


% Define psi angle range in degree
filament_psi_range = 45;



for t=1:189
    FilamentStruct2{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step2/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:189
    FilamentStruct3{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/

for k =1:189
    
    % load filament structure
    load(FilamentStruct2{k},'filament_struct_expanded');
    

    for i=1:size(filament_struct_expanded,2)
        filament_realindx = find(filament_struct_expanded(i).filament_processed_vec==1);
        filament_struct_expanded(i).filament_spin_cat(1:size(filament_struct_expanded(i).cor_filament_ext,1))= 0;
%         if size(filament_struct_expanded(i).cor_filament_ext,1)>2
        % extract every XY and fit
            filamentpoints = filament_struct_expanded(i).cor_filament_ext(:,1:3);
            b = [filamentpoints(end,:)-filamentpoints(1,:)]/norm([filamentpoints(end,:)-filamentpoints(1,:)]);
%             [b,stats] = robustfit(filamentpoints(:,1),filamentpoints(:,2));
%             b = polyfit(filamentpoints(:,1),filamentpoints(:,2),1);
        
%             scatter(filamentpoints(:,1),filamentpoints(:,2));
%             hold on
%             plot(filamentpoints(:,1),b(1)+b(2)*filamentpoints(:,1));
%             xlim([0 1024]);
%             ylim([0 1024]);
%             filament_struct_expanded(i).filament_robustfit_b = b;
%             filament_struct_expanded(i).filament_robustfit_stats = stats;

            Regresspsi = b;
            RegresspsiR = -b;
            filament_struct_expanded(i).slope = b;
            

            for j = 1:size(filament_realindx,1)
                seg_V = [];
                seg_V = filament_struct_expanded(i).rot_xyz_filament_ext(filament_realindx(j),:);
                if abs(acosd(dot(seg_V,Regresspsi)/(norm(seg_V)*norm(Regresspsi))))<=(filament_psi_range)
                    filament_struct_expanded(i).filament_spin_cat(filament_realindx(j)) = 1;
                elseif abs(acosd(dot(seg_V,RegresspsiR)/(norm(seg_V)*norm(RegresspsiR))))<=(filament_psi_range)
                    filament_struct_expanded(i).filament_spin_cat(filament_realindx(j)) = -1;
                else
                    filament_struct_expanded(i).filament_spin_cat(filament_realindx(j)) = 3;
                end
            end        
        


        
        
%         elseif size(filament_struct_expanded(i).cor_filament_ext,1) == 2
%             b = (filament_struct_expanded(i).cor_filament_ext(1,2)-filament_struct_expanded(i).cor_filament_ext(2,2))/...
%                 (filament_struct_expanded(i).cor_filament_ext(1,1)-filament_struct_expanded(i).cor_filament_ext(2,1));
%             Regresspsi = rad2deg(atan(b)) - 90;
%             RegresspsiR = wrapTo180(Regresspsi-180);
%             
%             for j = 1:2
%             if abs(angdiff(deg2rad(filament_struct_expanded(i).psi_filament_ext(j)),deg2rad(Regresspsi)))<=(pi/6)
%                 filament_struct_expanded(i).filament_spin_cat(j) = 1;
%             elseif abs(angdiff(deg2rad(filament_struct_expanded(i).psi_filament_ext(j)),deg2rad(RegresspsiR)))<=(pi/6)
%                 filament_struct_expanded(i).filament_spin_cat(j) = -1;
%             else
%                 filament_struct_expanded(i).filament_spin_cat(j) = 0;
%             end
%             end
%         end
            
        portA = size(find(filament_struct_expanded(i).filament_spin_cat==1),2); 
        portB = size(find(filament_struct_expanded(i).filament_spin_cat==-1),2);
        portmax = max([portA portB]);
        portconfi = portmax/size(filament_realindx,1);
%         if portconfi > 0.5
            if portA>portB
                filament_struct_expanded(i).filament_spin_indx = 1;
                filament_struct_expanded(i).spin_filament_conf = portconfi;
                filament_struct_expanded(i).mscS = portconfi*(1-portconfi)/size(filament_realindx,1);
                filament_struct_expanded(i).mscN = filament_struct_expanded(i).spin_filament_conf - 1.96*sqrt(filament_struct_expanded(i).mscS);
            elseif portA<portB
                filament_struct_expanded(i).filament_spin_indx = -1;
                filament_struct_expanded(i).spin_filament_conf = portconfi;
                filament_struct_expanded(i).mscS = portconfi*(1-portconfi)/size(filament_realindx,1);
                filament_struct_expanded(i).mscN = filament_struct_expanded(i).spin_filament_conf - 1.96*sqrt(filament_struct_expanded(i).mscS);
            else
                filament_struct_expanded(i).filament_spin_indx = [];
                filament_struct_expanded(i).spin_filament_conf = portconfi;
                filament_struct_expanded(i).mscS = portconfi*(1-portconfi)/size(filament_realindx,1);
                filament_struct_expanded(i).mscN = filament_struct_expanded(i).spin_filament_conf - 1.96*sqrt(filament_struct_expanded(i).mscS);
            end
%         end
            
    end
    

save(FilamentStruct3{k},'filament_struct_expanded');

disp(k)
end

clear

%% find regression line use fit_3d_data


% Define psi angle range in degree
filament_psi_range = 30;



for t=1:189
    FilamentStruct2{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step2/filaments_160_tomo_' num2str(t) '.mat'];
end

for t=1:189
    FilamentStruct3{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end

mkdir /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/

for k =1:189
    
    % load filament structure
    load(FilamentStruct2{k},'filament_struct_expanded');
    
    
    for i=1:size(filament_struct_expanded,2)
        
        % get processed data
        filament_realindx = find(filament_struct_expanded(i).filament_processed_vec==1);
        %initial filament_spin_cat
        filament_struct_expanded(i).filament_spin_cat(1:size(filament_struct_expanded(i).cor_filament_ext,1))= 0;
        % extract every XY and fit
        filamentpoints = filament_struct_expanded(i).cor_filament_ext(:,1:3);
        
        [coeffsxy, coeffsyx, coeffsxz, coeffszx, coeffsyz, coeffszy, b] = fit_3d_data(filamentpoints);
        
        filament_struct_expanded(i).slope = b;
        filament_struct_expanded(i).coeffsxy = coeffsxy;
        filament_struct_expanded(i).coeffsyx = coeffsyx;
        filament_struct_expanded(i).coeffsxz = coeffsxz;
        filament_struct_expanded(i).coeffszx = coeffszx;
        filament_struct_expanded(i).coeffsyz = coeffsyz;
        filament_struct_expanded(i).coeffszy = coeffszy;
        for j = 1:size(filament_realindx,1)
            seg_V = [];
            seg_V = filament_struct_expanded(i).rot_xyz_filament_ext(filament_realindx(j),:);
            Regresspsi = b(filament_realindx(j),:);
            RegresspsiR = -b(filament_realindx(j),:);
            if abs(acosd(dot(seg_V,Regresspsi)/(norm(seg_V)*norm(Regresspsi))))<=(filament_psi_range)
                filament_struct_expanded(i).filament_spin_cat(filament_realindx(j)) = 1;
            elseif abs(acosd(dot(seg_V,RegresspsiR)/(norm(seg_V)*norm(RegresspsiR))))<=(filament_psi_range)
                filament_struct_expanded(i).filament_spin_cat(filament_realindx(j)) = -1;
            elseif abs(acosd(dot(seg_V,Regresspsi)/(norm(seg_V)*norm(Regresspsi)))) <= abs(acosd(dot(seg_V,RegresspsiR)/(norm(seg_V)*norm(RegresspsiR))))
                filament_struct_expanded(i).filament_spin_cat(filament_realindx(j)) = 3;
            elseif abs(acosd(dot(seg_V,Regresspsi)/(norm(seg_V)*norm(Regresspsi)))) > abs(acosd(dot(seg_V,RegresspsiR)/(norm(seg_V)*norm(RegresspsiR))))
                filament_struct_expanded(i).filament_spin_cat(filament_realindx(j)) = -3;
            end
        end
        
        
        
        
        
        
        portA = size(find(filament_struct_expanded(i).filament_spin_cat==1),2);
        portB = size(find(filament_struct_expanded(i).filament_spin_cat==-1),2);
        portmax = max([portA portB]);
        portconfi = portmax/size(filament_realindx,1);
        %         if portconfi > 0.5
        if portA>portB
            filament_struct_expanded(i).filament_spin_indx = 1;
            filament_struct_expanded(i).spin_filament_conf = portconfi;
            filament_struct_expanded(i).mscS = portconfi*(1-portconfi)/size(filament_realindx,1);
            filament_struct_expanded(i).mscN = filament_struct_expanded(i).spin_filament_conf - 1.96*sqrt(filament_struct_expanded(i).mscS);
        elseif portA<portB
            filament_struct_expanded(i).filament_spin_indx = -1;
            filament_struct_expanded(i).spin_filament_conf = portconfi;
            filament_struct_expanded(i).mscS = portconfi*(1-portconfi)/size(filament_realindx,1);
            filament_struct_expanded(i).mscN = filament_struct_expanded(i).spin_filament_conf - 1.96*sqrt(filament_struct_expanded(i).mscS);
        else
            filament_struct_expanded(i).filament_spin_indx = [];
            filament_struct_expanded(i).spin_filament_conf = portconfi;
            filament_struct_expanded(i).mscS = portconfi*(1-portconfi)/size(filament_realindx,1);
            filament_struct_expanded(i).mscN = filament_struct_expanded(i).spin_filament_conf - 1.96*sqrt(filament_struct_expanded(i).mscS);
        end
        %         end
        
        
    end
    
    
    save(FilamentStruct3{k},'filament_struct_expanded');
    
    disp(k)
end

clear

%% add pseudo rot vector to the non-processed actin
for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end

for k =1:189
    
    % load filament structure
    load(FilamentStruct{k},'filament_struct_expanded');
    

    for i=1:size(filament_struct_expanded,2)
        coeffsxy = filament_struct_expanded(i).coeffsxy;
        coeffsyx = filament_struct_expanded(i).coeffsyx;
        coeffsxz = filament_struct_expanded(i).coeffsxz;
        coeffszx = filament_struct_expanded(i).coeffszx;
        coeffsyz = filament_struct_expanded(i).coeffsyz;
        coeffszy = filament_struct_expanded(i).coeffszy;
        for j = find(filament_struct_expanded(i).filament_processed_vec==0)'
            cor = [];
            cor = filament_struct_expanded(i).cor_filament_ext(j,:);
            b = [];
            b = psuedo_3d_data(cor,coeffsxy,coeffsyx,coeffsxz,coeffszx,coeffsyz,coeffszy);
            
            
            if filament_struct_expanded(i).filament_spin_indx == 1
                filament_struct_expanded(i).rot_xyz_filament_ext(j,:) = b;
                % find degree in y axis
                target_angle = filament_struct_expanded(i).rot_xyz_filament_ext(j,:);
                filament_struct_expanded(i).theta_filament_ext(j)=acosd(dot([0 1],[target_angle(1) target_angle(3)])/(norm([0 1])*norm([target_angle(1) target_angle(3)])));
                
                % find degree in psi after y rotation
                acosd(dot([0 1],[target_angle(1) target_angle(2)])/(norm([0 1])*norm([target_angle(1) target_angle(2)])));
                
                
            elseif filament_struct_expanded(i).filament_spin_indx == -1
                filament_struct_expanded(i).rot_xyz_filament_ext(j,:) = -b;
                % find degree in y axis
                target_angle = filament_struct_expanded(i).rot_xyz_filament_ext(j,:);
                filament_struct_expanded(i).rotate_yaxis(j) = acosd(dot([0 1],[target_angle(1) target_angle(3)])/(norm([0 1])*norm([target_angle(1) target_angle(3)])));
                
                % find degree in psi after y rotation
                filament_struct_expanded(i).rotate_zaxis(j) = acosd(dot([0.6104 0],[target_angle(1) target_angle(2)])/(norm([0 1])*norm([target_angle(1) target_angle(2)])));
            end

        end
        
        
        
        
    end
    
    
    save(FilamentStruct{k},'filament_struct_expanded');
    
    disp(k)
end

clear

%% plot check with regression line

for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end


mkdir /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step2_plot_2/

cd /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step2_plot_2/

for k = 1:189
    
    load(FilamentStruct{k},'filament_struct_expanded');
    
    figure(k)
    xlim([0 1024]);
    ylim([0 1024]);
    hold on
    for i =1:size(filament_struct_expanded,2)
        x = filament_struct_expanded(i).cor_filament_ext(:,1);
        y = filament_struct_expanded(i).cor_filament_ext(:,2);

        %         b = filament_struct_expanded(i).filament_robustfit_b;
        filamentpoints = filament_struct_expanded(i).cor_filament_ext(:,1:2);
        
        scatter(filamentpoints(:,1),filamentpoints(:,2),'.','b');
        
        %         plot(filamentpoints(:,1),b(1)+b(2)*filamentpoints(:,1),'r');
        
        if ~isempty( filament_struct_expanded(i).coeffsxy)
            plot(x, polyval(filament_struct_expanded(i).coeffsxy,x),'r');
            
        else
            plot(polyval(filament_struct_expanded(i).coeffsyx,y), y,'r');
            
        end
%         quiver(x,y,filament_struct_expanded(i).slope(:,1),filament_struct_expanded(i).slope(:,2))
        text(filamentpoints(size(filamentpoints,1),1),filamentpoints(size(filamentpoints,1),2),num2str(i),'Color','green','FontSize',8)
        
        
    end
    figname = ['tomo_' num2str(k)];
    set(gcf,'Position',[0 0 1024 1024],'Renderer','painters')
    
    saveas(gca, fullfile('/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step2_plot_2/',figname ),'png' )
    close(gcf)
    
    
    figure(k)
    xlim([0 1024]);
    ylim([0 750]);
    hold on
    for i =1:size(filament_struct_expanded,2)
        x = filament_struct_expanded(i).cor_filament_ext(:,1);
        y = filament_struct_expanded(i).cor_filament_ext(:,2);
        z = filament_struct_expanded(i).cor_filament_ext(:,3);
        %         b = filament_struct_expanded(i).filament_robustfit_b;
        filamentpoints = filament_struct_expanded(i).cor_filament_ext(:,[1 3]);
        
        scatter(x,z,'.','b');
        
        %         plot(filamentpoints(:,1),b(1)+b(2)*filamentpoints(:,1),'r');
        
        if ~isempty(filament_struct_expanded(i).coeffsxz)
            plot(x, polyval(filament_struct_expanded(i).coeffsxz,x),'r');
            
        elseif ~isempty( filament_struct_expanded(i).coeffszx)
            plot(polyval(filament_struct_expanded(i).coeffszx,z), z,'r');
            
        elseif ~isempty( filament_struct_expanded(i).coeffsyz)
            plot(x, polyval(filament_struct_expanded(i).coeffsyz,y),'r');
            
        elseif ~isempty( filament_struct_expanded(i).coeffszy)
            plot(x, z,'r');
            
        end
%         quiver(x,z,filament_struct_expanded(i).slope(:,1),filament_struct_expanded(i).slope(:,3))
        text(filamentpoints(size(filamentpoints,1),1),filamentpoints(size(filamentpoints,1),2),num2str(i),'Color','green','FontSize',8)
        
        
    end
    figname = ['tomo_side_xz' num2str(k)];
    set(gcf,'Position',[0 0 1024 1024],'Renderer','painters')
    
    saveas(gca, fullfile('/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step2_plot_2/',figname ),'png' )
    close(gcf)
    
    figure(k)
    xlim([0 1024]);
    ylim([0 750]);
    hold on
    for i =1:size(filament_struct_expanded,2)
        x = filament_struct_expanded(i).cor_filament_ext(:,1);
        y = filament_struct_expanded(i).cor_filament_ext(:,2);
        z = filament_struct_expanded(i).cor_filament_ext(:,3);
        %         b = filament_struct_expanded(i).filament_robustfit_b;
        filamentpoints = filament_struct_expanded(i).cor_filament_ext(:,[2 3]);
        
        scatter(y,z,'.','b');
        
        %         plot(filamentpoints(:,1),b(1)+b(2)*filamentpoints(:,1),'r');
        
        if ~isempty(filament_struct_expanded(i).coeffsxz)
            plot(y, polyval(filament_struct_expanded(i).coeffsxz,x),'r');
            
        elseif ~isempty( filament_struct_expanded(i).coeffszx)
            plot(y, z,'r');
            
        elseif ~isempty( filament_struct_expanded(i).coeffsyz)
            plot(y, polyval(filament_struct_expanded(i).coeffsyz,y),'r');
            
        elseif ~isempty( filament_struct_expanded(i).coeffszy)
            plot(polyval(filament_struct_expanded(i).coeffszy,z), z,'r');
            
        end
%         quiver(y,z,filament_struct_expanded(i).slope(:,2),filament_struct_expanded(i).slope(:,3))
        text(filamentpoints(size(filamentpoints,1),1),filamentpoints(size(filamentpoints,1),2),num2str(i),'Color','green','FontSize',8)
        
        
    end
    figname = ['tomo_side_yz' num2str(k)];
    set(gcf,'Position',[0 0 1024 1024],'Renderer','painters')
    
    saveas(gca, fullfile('/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step2_plot_2/',figname ),'png' )
    close(gcf)
end

clear
    
    
        
%% confidence score calculation

filament_com_conf_thres = 0.6;

for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end

for k = [1:28 30:110 113:189]
    load(FilamentStruct{k},'filament_struct_expanded');
    
    for i =1:size(filament_struct_expanded,2)
    
         % Calculate segment sensitivity confidence
         filament_struct_expanded(i).frac_of_pos_sense = size(find(filament_struct_expanded(i).segment_sensitivity_filament_ext==1),1)./size(find(filament_struct_expanded(i).filament_spin_cat~=0),2);
%          filament_struct_expanded(i).frac_of_pos_sense = size(find([filament_struct_expanded(i).segment_sensitivity_filament_ext==1]' & filament_struct_expanded(i).filament_spin_cat~=3 & filament_struct_expanded(i).filament_spin_cat~=0),2)./size(find(filament_struct_expanded(i).filament_spin_cat~=0),2);
         % Calculate combined confidence score / majority and sensitivity
         % filament_struct_ref(i).spin_filament_conf is [0.5,...,1.0]
         % filament_struct_ref(i).frac_of_pos_sense is  [0.0,...,1.0]
         filament_struct_expanded(i).filament_conf_com_score = filament_struct_expanded(i).spin_filament_conf .* filament_struct_expanded(i).frac_of_pos_sense;
    end
    
    
    % Select filaments based on combined confidence score
    filament_struct_ref_conf = [filament_struct_expanded.filament_conf_com_score];
    indx_bad = find(filament_struct_ref_conf < filament_com_conf_thres);
    disp(['Bundle: ' num2str(k) ' | Out of ' num2str(size(filament_struct_ref_conf,2)) ' filaments in this bundle ' num2str(size(indx_bad,2)) ' will be excluded | Number of filaments: ' num2str(size(filament_struct_ref_conf,2)-size(indx_bad,2))]);
%     filament_struct_expanded(indx_bad) = [];
    clear filament_struct_ref_conf;
    clear indx_bad;
save(FilamentStruct{k},'filament_struct_expanded');

disp(k)    
end

clear


%% plot ccs score 
for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end


mkdir /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/multiplot/

cd /home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/multiplot/



score = [];

for k = [1:28 30:110 113:189]
    load(FilamentStruct{k},'filament_struct_expanded');

    for i = 1:size(filament_struct_expanded,2)
    if ~isempty(filament_struct_expanded(i).filament_spin_indx)
        score = [score filament_struct_expanded(i).filament_conf_com_score];

    end
    end
    
end
histogram(score,40,'EdgeColor','none')
good = find(score>=0.6);
goodR = size(good,2)/size(score,2);

set(gcf,'Position',[0 0 1024 1024],'Renderer','painters')
box off


for k = 1:189
    subplot(8,8,k)
    load(FilamentStruct{k},'filament_struct_expanded');
    
    score = [filament_struct_expanded.filament_conf_com_score];
    good = find(score>0.6);
    goodR = size(good,2)/size(score,2);
    histogram(score,20,'EdgeColor','none')
    title([sprintf('tomo_ %d',k),' : ',num2str(goodR)]);
    xticks([0 0.6 1])

    
    
end
set(gcf,'Position',[0 0 1024 1024],'Renderer','painters')

saveas(gcf,'ccs_score','tiffn')
        
f = figure(2);

for k = 1:58
    subplot(8,8,k)
    load(FilamentStruct{k},'filament_struct_expanded');
    
    score = [filament_struct_expanded.mscN];
    good = find(score>0.5);
    goodR = size(good,2)/size(score,2);
    histogram(score,20,'EdgeColor','none')
    title([sprintf('tomo_ %d',k),' : ',num2str(goodR)]);
    xticks([0 0.5 1])

    
    
end
set(gcf,'Position',[0 0 1024 1024],'Renderer','painters')

saveas(gcf,'msc_score','tiffn')      


f = figure(3);


for k = 1:58
    subplot(8,8,k)
    hold on
    load(FilamentStruct{k},'filament_struct_expanded');
    score = [];
    length = [];
    for i =1:size(filament_struct_expanded,2)
    
        score = [score filament_struct_expanded(i).mscN];
        length = [length size(find([filament_struct_expanded(i).filament_processed_vec]~=0),1)];
    
    end
    scatter(length,score,7,'filled')
        

     title(sprintf('tomo_ %d',k))
     yticks([0 0.5 1])
    
    
end        
  
set(gcf,'Position',[0 0 1024 1024],'Renderer','painters')

saveas(gcf,'msc2length_score','tiffn')


% figure(1)
% hold on
% 
% figure(2)
% hold on
% for k = 2: 58
%     load(FilamentStruct{k},'filament_struct_expanded');
%     scoreM = [];
%     scoreC = [];
%     for i =1:size(filament_struct_expanded,2)
%         scoreM = [scoreM filament_struct_expanded(i).mscN];
%         scoreC = [scoreC filament_struct_expanded(i).filament_conf_com_score];
%     end
%     figure(1)
%     cdfplot(scoreM)
%     figure(2)
%     cdfplot(scoreC)
% end
    
        
       
close all
clear       
                
            
%% plot after ccs filtered

%withouht color
for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end

for k =[1:28 30:110 113:189]
    
    figure(k)
    hold on
    xlim([0 1024])
    ylim([0 1024])
    zlim([0 1024])
    rrr = [];
    kkk = [];
    load(FilamentStruct{k},'filament_struct_expanded');
    for i = 1:size(filament_struct_expanded,2)
        
        if filament_struct_expanded(i).filament_conf_com_score >= 0.6
            
            filament_process = [];
            filament_process = filament_struct_expanded(i).filament_processed_vec;
            rrr = [rrr;filament_struct_expanded(i).cor_filament_ext(filament_process==1,:)];
            kkk = [kkk;filament_struct_expanded(i).cor_filament_ext(filament_process==0,:)];
            
            plot3(filament_struct_expanded(i).cor_filament_ext(:,1)...
                ,filament_struct_expanded(i).cor_filament_ext(:,2)...
                ,filament_struct_expanded(i).cor_filament_ext(:,3),'k')
            
            text(filament_struct_expanded(i).cor_filament_ext(1,1)...
                ,filament_struct_expanded(i).cor_filament_ext(1,2)...
                ,filament_struct_expanded(i).cor_filament_ext(1,3),num2str(i),'Color','green','FontSize',8)
            
            
            
            
        end
        
        
        
        
        
        
    end
    
    scatter3(rrr(:,1),rrr(:,2),rrr(:,3),'.','r')
%     scatter3(kkk(:,1),kkk(:,2),kkk(:,3),'*','k')
    figname = ['tomo_new' num2str(k)];
    set(gcf,'Position',[1 1 1024 1024],'Renderer','painters')
    
    saveas(gca, fullfile('/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3_plot_2/',figname ),'png' )
    close(gcf)
    
    
    
end


clear

            
% with color
for t=1:189
    FilamentStruct{t} = ['/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3/filaments_160_tomo_' num2str(t) '.mat'];
end

for k =1:189
    
    figure(k)
    hold on
    xlim([0 1024])
    ylim([0 1024])
    zlim([0 1024])
    turbo_wrap = turbo(361);
    colormap(turbo_wrap)
    
    load(FilamentStruct{k},'filament_struct_expanded');
    for i = 1:size(filament_struct_expanded,2)
        
        if filament_struct_expanded(i).filament_conf_com_score >= 0.6
            
            filament_process = [];
            filament_process = filament_struct_expanded(i).filament_processed_vec;
            psi = round(filament_struct_expanded(i).psi_filament_ext(filament_process==1)+180);
            [idx,C] = kmeans(psi,2,'Replicates',20);
            if size(find(idx==1),1) >= size(find(idx==2),1)
                psi(idx==2) = wrapTo360(psi(idx==2)-180);
            else
                psi(idx==1) = wrapTo360(psi(idx==1)-180);
            end
            
            scatter3(filament_struct_expanded(i).cor_filament_ext(filament_process==1,1),...
                filament_struct_expanded(i).cor_filament_ext(filament_process==1,2),...
                filament_struct_expanded(i).cor_filament_ext(filament_process==1,3),15,turbo_wrap(psi+1,:),'filled')
            
            plot3(filament_struct_expanded(i).cor_filament_ext(:,1)...
                ,filament_struct_expanded(i).cor_filament_ext(:,2)...
                ,filament_struct_expanded(i).cor_filament_ext(:,3),'k')
            
            %             text(filament_struct_expanded(i).cor_filament_ext(1,1)...
            %                ,filament_struct_expanded(i).cor_filament_ext(1,2)...
            %                ,filament_struct_expanded(i).cor_filament_ext(1,3),num2str(i),'Color','green','FontSize',8)
            
            
            
            
        end
        
        
        
        
        
        
    end
    colorbar;
    colorbar('Ticks',[0:0.25:1],'TickLabels',{'0','90','180','270','360'})
    figname = ['tomo_' num2str(k)];
    set(gcf,'Position',[1 1 1024 1024],'Renderer','painters')
    
    saveas(gca, fullfile('/home/Medalia/Projects7/WChung/actin_adhesion/mapping3d/filaments_step3_plot_2_colored/',figname ),'png' )
    close(gcf)
    
    
    
end


clear
                        
  
            

        
        
