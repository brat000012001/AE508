function plots(t,X,r0,rf,mu)

    function plot_control(t,X)
        u = X(:,11:13)./vecnorm(X(:,11:13),2,2);
    
        figure;
        title('Control');
        subplot(3,1,1);
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


    function plot_trajectory(t, X, r0, rf, mu)
        figure;  
        plot3(X(:,1),X(:,2),X(:,3));
        hold on;
        plot3(r0(1),r0(2),r0(3),'*b');
        plot3(rf(1),rf(2),rf(3),'*g');
        %plot_asteroid_trajectory(t(1), t(end), r0, v0, mu);
        hold off;
    end


    function plot_costates(t, X)
        figure;
        grid on; hold on;
    
        title('Costate \lambda_r');
        subplot(3,1,1);
        plot(t, X(:,8),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_r_x');
    
        subplot(3,1,2);
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
        grid on; hold on;
        title ('Costate \lambda_v');
        subplot(3,1,1);
        plot(t, X(:,11),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_v_x');
    
        subplot(3,1,2);
        plot(t, X(:,12),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_v_y');
    
        subplot(3,1,3);
        plot(t, X(:,13),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_v_z');
        hold off;
    
        % Costate lambda_m
        figure;
        grid on; hold on;
        title('\lambda_m');
        plot(t, X(:,14),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('\lambda_m');
        hold off;
    end


    function plot_states(t,X)
        figure;
        grid on; hold on;
    
        title("Position")
        subplot(3,1,1);
        plot(t, X(:,1),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('x (km)');
    
        subplot(3,1,2);
        plot(t, X(:,2),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('y (km)');
        
        subplot(3,1,3);
        plot(t, X(:,3),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('z (km)');
    
        hold off;
    
        % Velocity
        figure;
        grid on; hold on;
        
        title("Velocity")
        subplot(3,1,1);
        plot(t, X(:,4),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('v_x (km/s)');
    
        subplot(3,1,2);
        plot(t, X(:,5),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('v_y (km/s)');
        
        subplot(3,1,3);
        plot(t, X(:,6),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('v_z (km/s)');
        
        hold off;
    
        % Mass 
        figure;
        grid on; hold on;
    
        title("Mass")
        plot(t, X(:,7),'LineWidth',1.5);
        xlabel('t (sec)');
        ylabel('m (kg)');
        
        hold off;
    end

    %{
    function plot_asteroid_trajectory(t0,tf,r0,v0,mu)
    
        function XDot = eom_aroid(t,X,mu)
            r = X(1:3);
            v = X(4:6);
            rdot = v;
            vdot = -mu/norm(r)^3*r;
            XDot = [rdot; vdot];
        end
    
        opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
        [t, X] = ode45(@eom_aroid, [t0 tf], [r0; v0],opts_ode,mu);
        plot3(X(:,1),X(:,2),X(:,3),'LineWidth',1.5);
    end
    %}

    plot_trajectory(t,X,r0,rf,mu);
    plot_states(t,X);
    plot_costates(t,X);
    plot_control(t,X);

end