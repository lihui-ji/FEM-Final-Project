close all
clear all
clc
set(groot,'defaultFigureVisible','off')
cd 'C:\DUKE\courses\FEM\project'
%% Convergence by mesh refinement
meshsizes=[2,5,10,25];%mesh size
convergence_mesh=[];
resonance_case=1 %set 0 for normal case
if resonance_case
    n_resonance=1;%resonance eigenmode in x
    m_resonance=1;%resonance eigenmode in y
    omega=sqrt(n_resonance*pi^2+m_resonance*pi^2); %frequency for eigen value lambda_nm;
else
    omega=1;
end
for ih=1:4
    clearvars -except meshsizes convergence_mesh ih omega n_resonance m_resonance convergence_time it time_steps
    m=meshsizes(ih);
    t=5;
    dt=1;%time step
    Ndt=t/dt;%maximum iteration time
    dt_save=2;%time step save interval
    savecount=0;
    t_iteration=0;
    FEM_final_project_maincode;
    % analytical solution
    syms funx funy funt
    for m_mode=1:20
        for n_mode=1:20
            resonance=(n_mode==n_resonance)*(m_mode==m_resonance);
            phi(m_mode,n_mode)=sin(n_mode*pi*funx).*sin(m_mode*pi*funy);
            lambda(m_mode,n_mode)=(n_mode*pi)^2+(m_mode*pi)^2;
            gamma(m_mode,n_mode)=4*(cos(n_mode*pi)-1)*(cos(m_mode*pi)-1)/m_mode/n_mode/pi/pi;
            if resonance
                A(m_mode,n_mode)=gamma(m_mode,n_mode)*funt*sin(omega*funt)/2/omega;
            else
                A(m_mode,n_mode)=gamma(m_mode,n_mode)*cos(omega*funt)/(lambda(m_mode,n_mode)-omega^2);
            end
            u_mode(m_mode,n_mode)=A(m_mode,n_mode).*phi(m_mode,n_mode);        
        end
    end
    u=matlabFunction(sum(sum(u_mode)),'vars',[funx funy funt]);

    for i=1:num_n
        u_analytical(i)=u(node_coord(1,i),node_coord(2,i),t);
    end
    area=0;
    for i = 1:length(num_e)
        for j = 1:3
            node_index=ele_conn(j,i);
            area=area+(uh_new(node_index)-u_analytical(node_index))^2*Se/3;
        end
    end
    convergenceh=sqrt(area);
    convergence_mesh(ih)=convergenceh;
end
%% Convergence by time refinement
time_steps=[1 0.5 0.25 0.125];
convergence_time=[];
resonance_case=1; %set 0 for normal case
if resonance_case
    n_resonance=1;%resonance eigenmode in x
    m_resonance=1;%resonance eigenmode in y
    omega=sqrt(n_resonance*pi^2+m_resonance*pi^2); %frequency for eigen value lambda_nm;
else
    omega=1;
end
for it=1:4
    clearvars -except meshsizes convergence_mesh ih omega n_resonance m_resonance convergence_time it time_steps
    m=25;
    t=5;
    dt=time_steps(it);%time step
    Ndt=t/dt;%maximum iteration time
    dt_save=2;%time step save interval
    savecount=0;
    t_iteration=0;
    FEM_final_project_maincode;
    % analytical solution
    syms funx funy funt
    for m_mode=1:20
        for n_mode=1:20
            resonance=(n_mode==n_resonance)*(m_mode==m_resonance);
            phi(m_mode,n_mode)=sin(n_mode*pi*funx).*sin(m_mode*pi*funy);
            lambda(m_mode,n_mode)=(n_mode*pi)^2+(m_mode*pi)^2;
            gamma(m_mode,n_mode)=4*(cos(n_mode*pi)-1)*(cos(m_mode*pi)-1)/m_mode/n_mode/pi/pi;
            if resonance
                A(m_mode,n_mode)=gamma(m_mode,n_mode)*funt*sin(omega*funt)/2/omega;
            else
                A(m_mode,n_mode)=gamma(m_mode,n_mode)*cos(omega*funt)/(lambda(m_mode,n_mode)-omega^2);
            end
            u_mode(m_mode,n_mode)=A(m_mode,n_mode).*phi(m_mode,n_mode);        
        end
    end
    u=matlabFunction(sum(sum(u_mode)),'vars',[funx funy funt]);

    for i=1:num_n
        u_analytical(i)=u(node_coord(1,i),node_coord(2,i),t);
    end
    area=0;
    for i = 1:length(num_e)
        for j = 1:3
            node_index=ele_conn(j,i);
            area=area+(uh_new(node_index)-u_analytical(node_index))^2*Se/3;
        end
    end
    convergencet=sqrt(area);
    convergence_time(it)=convergencet;
end
%% visualize results mesh refinement
figure;
loglog(1/4./meshsizes./meshsizes,convergence_mesh,'o')
p=polyfit(log(1/4./meshsizes./meshsizes),log(convergence_mesh),1);
hold on
xi=logspace(-5,0,5);
loglog(xi,exp(polyval(p,log(xi))))
title({'Convergence Result for h refinement',['slope is ',num2str(p(1))]})
xlabel('h')
ylabel('error')
grid on
saveas(gcf,['Convergence Result for h refinement at t=',num2str(t),'s ','resonance','.png'])
%% visualize results time refinement
figure;
plot(time_steps,convergence_time,'o')
p=polyfit(time_steps,convergence_time,1);
slope=polyfit(log(time_steps),log(convergence_time),1);
hold on
xi=linspace(0,1);
plot(xi,polyval(p,xi))
title({'Convergence Result for time refinement','first order accuracy in time'})
xlabel('time step/s')
ylabel('error')
grid on
saveas(gcf,['Convergence Result for time refinement at m=',num2str(m),' ','resonance','.png'])