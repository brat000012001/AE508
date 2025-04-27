function S = plots(mu)
    addpath(".");

    function plot_control(t,X)
        u = X(:,11:13)./vecnorm(X(:,11:13),2,2);
    
        figure; hold on; grid on;
        sgtitle('Control u');

        subplot(3,1,1); hold on;
        plot(t,u(:,1),'LineWidth', 1.5);
        xlabel('t (sec)');
        ylabel('u_x');
    
        subplot(3,1,2);
        plot(t,u(:,2),'LineWidth', 1.5);
        xlabel('t (sec)');
        ylabel('u_y');
        
        subplot(3,1,3);
        plot(t,u(:,3),'LineWidth', 1.5);
        xlabel('t (sec)');
        ylabel('u_z');
    end


    function plot_trajectory(t, X, r0, vgeo, rf, vf)
        figure;
        plot3(X(:,1),X(:,2),X(:,3),'b-','LineWidth',1.5);
        hold on;
        plot3(0,0,0,'b*','LineWidth',3);
        plot3(r0(1),r0(2),r0(3),'g*','LineWidth',3);
        plot3(rf(1),rf(2),rf(3),'r*','LineWidth',3);

        % Compute and plot a solution to Lambert's problem
        [V1,V2] = lambert(r0,rf,t(end),mu,'pro');
        opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
        [Tlambert, Xlambert] = ode45(@two_body, [t(1) t(end)], ...
            [r0; V1], opts_ode, mu);
        [~, Xast] = ode45(@two_body, [0 86400], ...
            [rf; vf], opts_ode, mu);
        [~, Xgeo] = ode45(@two_body, [0 86400], ...
            [r0; vgeo], opts_ode, mu);
        
        plot3(Xlambert(:,1),Xlambert(:,2),Xlambert(:,3),'m--','LineWidth',1.5);
        plot3(Xast(:,1),Xast(:,2),Xast(:,3),'r--');
        plot3(Xgeo(:,1),Xgeo(:,2),Xgeo(:,3),'b--');
        legend('Optimal Trajectory', ...
            'Earth Center',...
            'SC @ 0',...
            'AST @ t_f',...
            'Lambert Solution',...
            'Asteroid Trajectory @ t_f',...
            'Departure Orbit @ t0');
        title(sprintf('Optimal (v=%.3f) vs Lambert (v2 = %.3f)',norm(X(end,4:6)),norm(V2)));
        xlabel('x (km)');
        ylabel('y (km)');
        zlabel('z (km)');
        hold off;
    end


    function plot_costates(t, X)
        figure;
        grid on; hold on;
        sgtitle('\lambda_r');

        subplot(3,1,1); hold on;
        plot(t, X(:,8),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_r_x');
    
        subplot(3,1,2); hold on;
        plot(t, X(:,9),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_r_y');
    
        subplot(3,1,3);
        plot(t, X(:,10),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_r_z');
        hold off;
        
        % lambda_v
        figure;
        sgtitle('\lambda_v');
        grid on; hold on;
        subplot(3,1,1);
        plot(t, X(:,11),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_v_x');
    
        subplot(3,1,2); hold on;
        plot(t, X(:,12),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_v_y');
    
        subplot(3,1,3); hold on;
        plot(t, X(:,13),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_v_z');
    
        % Costate lambda_m
        figure;
        grid on; hold on;
        sgtitle('\lambda_m');
        plot(t, X(:,14),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_m');
    end


    function plot_states(t,X)
        figure;
        grid on; hold on;
        sgtitle('Position km');

        subplot(3,1,1); hold on;
        plot(t, X(:,1),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('x (km)');
    
        subplot(3,1,2); hold on;
        plot(t, X(:,2),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('y (km)');
        
        subplot(3,1,3); hold on;
        plot(t, X(:,3),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('z (km)');
    
        % Velocity
        figure;
        grid on; hold on;
        sgtitle('Velocity km/s');
        
        subplot(3,1,1); hold on;
        plot(t, X(:,4),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('v_x (km/s)');
    
        subplot(3,1,2); hold on;
        plot(t, X(:,5),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('v_y (km/s)');
        
        subplot(3,1,3); hold on;
        plot(t, X(:,6),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('v_z (km/s)');
    
        % Mass 
        figure;
        grid on; hold on;
        sgtitle('Mass kg');
    
        plot(t, X(:,7),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('m (kg)');

    end

    S.trajectory = @(t,X,r0,vgeo,rf,vf) plot_trajectory(t,X,r0,vgeo,rf,vf);
    S.states = @(t,X) plot_states(t,X);
    S.costates = @(t,X) plot_costates(t,X);
    S.control = @(t,X) plot_control(t,X);

end