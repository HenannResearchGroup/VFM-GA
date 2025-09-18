function [lambda_out] = S_sub_evaluate(thetas)
    fileid = fopen('record_expdata_name.txt');
    expdata_name = textscan(fileid,'%s','delimiter','\n');
    expdata_name = expdata_name{1}{1};
    load(expdata_name);
    
    record_input_setup = fopen('record_input_data_setup.txt','r');
    info_temp = textscan(record_input_setup,'%d %d %d\n');
    data_mod = info_temp{1};
    
    
    

    %% input parameters
    min_e11 = min(compression_ax_strain);
    if data_mod == 1 || data_mod == 3
        max_e11 = max(compression_ax_strain);
    end
    time1 = 1500;
    time2 = 1000;
    step = 500;
    orange = 	'#EDB120';
    
    
    fig_prop = {'k-','linewidth',1};
    figure(1)
    hold on
    plot(compression_ax_strain,compression_ax_stress,'k-','LineWidth',3)
    hold on
    if data_mod == 1 || data_mod == 3
        plot(tension_ax_strain,tension_ax_stress,'k-','LineWidth',3)
        xlim([min(compression_ax_strain)-0.05 max(tension_ax_strain)+0.05])
        ylim([min(compression_ax_stress)*1.2 max(tension_ax_stress)*1.2])
    else
        xlim([min(compression_ax_strain)-0.05 0])
        ylim([min(compression_ax_stress)*1.2 0])
    end
    hold on 

    set(gca,'FontSize',20)
    xlabel('Axial engineering strain','fontsize',19)
    ylabel('Axial engineering stress [{\itPa}]','fontsize',19)
    
    if data_mod == 2 || data_mod == 1 
        figure(2)
        hold on 
        plot(compression_ax_strain,compression_lat_strain,'k-','LineWidth',3)
        hold on 
        if data_mod == 1 
            plot(tension_ax_strain,tension_lat_strain,'k-','LineWidth',3)
            xlim([min(compression_ax_strain)-0.05 max(tension_ax_strain)+0.05])
            ylim([min(tension_lat_strain)-0.1 max(compression_lat_strain)+0.1])
        else
            xlim([min(compression_ax_strain)-0.05 0])
            ylim([0 max(compression_lat_strain)+0.1])
        end
        hold on
        
        
        set(gca,'FontSize',20)
        xlabel('Axial engineering strain','fontsize',19)
        ylabel('Lateral engineering strain','fontsize',19)
    end
    %% construct whole thetas with full parameters, but actually not
    %needed
    if length(thetas) ~= 14
        display('The material parameter input is wrong');
        display('Only plotting exp data now, go back and replugin the correct mat parameters');
    else
        theta(:) = thetas(:);

        PROPS155 = [theta,...
            0,0,0,0,0,0,...
            0,0,0,0,0,0,...
            0,...
            0,0,...%Gneq1,t1
            0,0.,...%Gneq2,t2
            0,0,...%Gneq3,t3
            0,0,...
            0, 0, 0, 0];
        PROPS = PROPS155;
        types = [1, 1, 1, 1]; %flow rule type for each branch, 2 is shear thining, 1 is newtonian 
        branches = [0, 0, 0, 0];
        
        %% compute
        step2 = step-1;
        increment_t = time1/step2;
        strain_t = min_e11/step2*2;
        e11_mean = [0:strain_t:min_e11,min_e11:-strain_t:0];

        time_mean = 0:increment_t:time1;
        
        dt = time1/step;
        time = 0:dt:time1;

        lambda1 = spline(time_mean, e11_mean, time);
        lambda1 = lambda1+1;
        lambda_out = lambda1;
        [stress1, lambda2, K1, ~]  = S_subsub_solve_explicit(lambda1, dt, PROPS, types, branches);

        if data_mod == 1 || data_mod == 2
            figure(2)
            hold on
            plot(lambda1-1, lambda2-1,'b--','linewidth',3)
        end
        figure(1)
        hold on
        plot(lambda1-1, stress1,'b--','linewidth',3)
    if data_mod == 1 || data_mod == 3 

          theta(:) = thetas(:);
            PROPS155 =  [theta,...
            0,0,0,0,0,0,...
            0,0,0,0,0,0,...
            0,...
            0,0,...%Gneq1,t1
            0,0.,...%Gneq2,t2
            0,0,...%Gneq3,t3
            0,0,...
            0, 0, 0, 0];
        PROPS = PROPS155;
        types = [1, 1, 1, 1]; %flow rule type for each branch, 2 is shear thining, 1 is newtonian 
        branches = [0, 0, 0, 0];
        %% start computing
        step2 = step-1;
        increment_t = time2/step2;
        strain_t = max_e11/step2*2;
        e11_mean = [0:strain_t:max_e11,max_e11:-strain_t:0];
        time_mean = 0:increment_t:time2;
        dt = time2/step;
        time = 0:dt:time2;
        lambda1 = spline(time_mean, e11_mean, time);
        lambda1 = lambda1+1;
        [stress1, lambda2,  K1,~ ]  = S_subsub_solve_explicit(lambda1, dt, PROPS, types, branches);
        if data_mod == 1 
            figure(2)
            hold on
            plot(lambda1-1, lambda2-1, 'b--', 'linewidth',3)
        end
        figure(1)
        hold on
        plot(lambda1-1, stress1,'b--', 'linewidth',3)
    end
    end


end
