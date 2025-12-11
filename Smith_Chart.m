clc
clear
close all

%We will start by defining the data given on the presentation for our case,
%HP turbine (SLIDE 7), plus the initial values needed for the viscosity
%calculations (SLIDE 23). All units are in IS.

m_dot=150;
N=14000*pi/30;
T01=2000;
p01=45e5;
gamma=1.3;
R_gas=287;
varT0_T01=0.17;
rmax=0.40;
R=0.6;
cp=gamma*R_gas/(gamma-1);
M_max=0.75;
eps_max=95*pi/180;
VR_min=1.5; %c2/c1=> stator ; w3/w2=> rotor
zweifel=0.8;
t_over_c=0.2;
H_over_b=4;
viscosity0=1.7894e-5;
T0=273.11;
S=110.56;


%We will now make an initial choice for the mean radius and calculate U,
%speed of sound and Mach for our choice

rm_init=0.36;

U_init=N*rm_init;
a01=sqrt(gamma*R_gas*T01);
Mr=U_init/a01;

%We will start by making a mesh for both, the loading and the flow
%coefficients and getting the values for the velocity triangles for our
%choice. 

[phi, psi]=meshgrid(linspace(0.2, 1.5, 200), linspace(0.2, 3, 200));
[alpha_1,alpha_2,beta_2,beta_3,epsilon_s,epsilon_r,c1_U,c1,c2_U,c2,w2_U,w2,w3_U,w3]=turbine_triangle(phi,psi,R,U_init);

%In order to fully plot the Smith chart we need to calculate the efficiency
%for the different loading and flow coefficients. The formula, given on
%SLIDE 13, will require the calculation of several temperature ratios and
%the loss coefficients

T2_T01=1-(gamma-1)/2* Mr^2.*(c2_U.^2);
a2=sqrt(gamma*R_gas.*(T2_T01*T01));
M2=(c2_U*U_init)./a2;

T03_T01=1-(gamma-1)*Mr^2.*psi;

T3_T01=T03_T01-(gamma-1)/2*Mr^2.*(c1_U.^2);
a3_rel=sqrt(gamma*R_gas.*(T3_T01.*T01));
M3_rel=(w3_U.*U_init)./a3_rel;

%The loss coefficients (SLIDE 15) will need to be corrected for the aspect
% ratio, as for our case it is larger than 3 (4, to be exact), and we will
% need to check if it should be corrected for the Reynolds

loss_coeff_r=struct("star",0, "Hb_corr", 0, "Re_corr", 0, "star_opt", 0, ...
    "Hb_corr_opt", 0, "Re_corr_opt", 0);
loss_coeff_s=struct("star",0, "Hb_corr", 0, "Re_corr", 0, "star_opt", 0, ...
    "Hb_corr_opt", 0, "Re_corr_opt", 0);

loss_coeff_s.star=0.04+0.06*(epsilon_s*180/(pi*100)).^2;
loss_coeff_r.star=0.04+0.06*(epsilon_r*180/(pi*100)).^2;

loss_coeff_s.Hb_corr=(1+loss_coeff_s.star)*(0.993+0.021/H_over_b)-1;
loss_coeff_r.Hb_corr=(1+loss_coeff_r.star)*(0.975+0.075/H_over_b)-1;

%To calculate the Re, we need to get the values of rho2, rho3, and the
%different Dh, in addition to the values of mu2 and mu3 (SLIDE 23) as the
%viscosity depends on the static temperature 

p2=p01.*T2_T01.^(gamma/(gamma-1));   % isentropic
rho2=p2./(R_gas.*T01.*T2_T01);

p3=p01.*T3_T01.^(gamma/(gamma-1));   % isentropic
rho3=p3./(R_gas.*T01.*T3_T01);

T1=T01-c1.^2/(2*cp);
p1=p01./(T01./T1).^(gamma/(gamma-1));
rho1=p1./(R_gas.*T1);

Area=m_dot./(rho2*phi.*U_init);
H=Area./(2*pi*rm_init);

s_over_b_s=0.4 ./(cos(alpha_2).^2.*(tan(alpha_2)+tan(alpha_1)));
s_over_b_r=0.4 ./(cos(beta_3).^2.*(tan(beta_2)+tan(beta_3)));

b=H/H_over_b;
s_s=s_over_b_s.*b;
s_r=s_over_b_r.*b;

Dh_s=(2*s_s.*H.*cos(alpha_2))./(s_s.*cos(alpha_2)+H); 
Dh_r=(2*s_r.*H.*cos(beta_3))./(s_r.*cos(beta_3)+H);

viscosity2=viscosity0.*((T2_T01*T01)./T0).^(3/2).* ((T0+S)./(T2_T01*T01+S));
viscosity3=viscosity0.*((T3_T01*T01)./T0).^(3/2).* ((T0+S)./(T3_T01*T01+S));

Re_s=rho2.*c2.*Dh_s./viscosity2;
Re_r=rho3.*w3.*Dh_r./viscosity3;

%As the Re is different of 10^5, we need to correct for the Re
loss_coeff_s.Re_corr=(10^5./Re_s).^0.25 .*loss_coeff_s.Hb_corr;
loss_coeff_r.Re_corr=(10^5./Re_r).^0.25 .*loss_coeff_r.Hb_corr;

%As the formula for the density has several terms, we will split it to make
%the calculations easier and reduce the risk of mistakes

bracket_s=T2_T01.^-1 .*(loss_coeff_s.Re_corr./(1-loss_coeff_s.Re_corr)).*c2_U.^2;
bracket_r=T3_T01.^-1 .*(loss_coeff_r.Re_corr./(1-loss_coeff_r.Re_corr)).*w3_U.^2;
denominator=1+(1./(2.*psi)).*T03_T01.*(bracket_s+bracket_r);
efficiency=1./denominator;

%We will now plot the Smith chart and draw the boundary lines for our case

contour(phi,psi,efficiency)   %ADD COLORS TO CONTOUR
hold on
xlabel("Flow factor, \phi")
ylabel("Loading factor, \psi")

%So that we can extract the best value for the flow and loading
%coefficients we first need to plot the limits imposed in the guide (SLIDE
%7)

VR_s=c2./c1;
VR_r=w3./w2;

contour(phi, psi, M2, [M_max M_max],"Color","red","Linewidth",1.2);
contour(phi, psi, M3_rel, [M_max M_max], "Color","#e69e19","Linewidth",1.2);
contour(phi, psi, VR_s, [VR_min VR_min], "Color","green","Linewidth",1.2); 
contour(phi, psi, VR_r, [VR_min VR_min],"Color","#0f430f","Linewidth",1.2); 
contour(phi, psi, epsilon_s, [eps_max eps_max], "Color","#808080","Linewidth",1.2);
contour(phi, psi, epsilon_r, [eps_max eps_max], "Color","black","Linewidth",1.2);

%Now that we have the contour we will highlight the zone in which all the
%conditions are fullfiled so that we can make a choice for the most
%appropiate loading and flow coefficients. (Made with ChatGPT)

valid=(M2 < M_max) & (M3_rel < M_max) & (VR_s > VR_min) & (VR_r > VR_min) & ...
        (epsilon_s < eps_max) & (epsilon_r < eps_max);
contourf(phi, psi, valid, [0.5 1], "FaceAlpha", 0.3, "LineColor", "none", "FaceColor", "cyan");


%Choose the phi and psi

phi_opt=0.7;
psi_opt=1.5;
plot(phi_opt, psi_opt, "x");
legend("\phi and \psi contours","M2", "M3_{rel}","VR_s","VR_r","\epsilon_s",...
    "\epsilon_r", "Design area","Point chosen", "Location","best")
hold on

stages=(cp*T01*varT0_T01)/(psi_opt*U_init^2);
stages=round(stages);

% Calculate the velocity triangles for the optimum values of psi and phi

[alpha_1_opt,alpha_2_opt,beta_2_opt,beta_3_opt,epsilon_s_opt,epsilon_r_opt,...
    c1_U_opt,c1_opt,c2_U_opt,c2_opt,w2_U_opt,w2_opt,w3_U_opt,w3_opt...
    ]=turbine_triangle(phi_opt,psi_opt,R,U_init);

var_h0=psi_opt*U_init^2;
var_T0=var_h0/cp;  %MIRAR SI HAY QUE DIVIDIR ENTRE LOS STAGES

loss_coeff_s.star_opt=0.04+0.06*(epsilon_s_opt*180/(pi*100))^2;
loss_coeff_r.star_opt=0.04+0.06*(epsilon_r_opt*180/(pi*100))^2;

loss_coeff_s.Hb_corr_opt=(1+loss_coeff_s.star_opt)*(0.993+0.021/H_over_b)-1;
loss_coeff_r.Hb_corr_opt=(1+loss_coeff_r.star_opt)*(0.975+0.075/H_over_b)-1;

for i=1:stages

    if i==1
        T01_opt(i)=T01;
        p01_opt(i)=p01;
    else
        T01_opt(i)=T03_opt(i-1); 
        p01_opt(i)=p03_opt(i-1);
    end

    T1_opt(i)=T01_opt(i)-c1_opt^2/(2*cp);
    p1_opt(i)=(T1_opt(i)/T01_opt(i))^(gamma/(gamma-1))*p01_opt(i); %MIRAR T1S
    Mr_opt(i)=U_init/sqrt(gamma*R_gas*T01_opt(i));
    rho1_opt(i)=p1_opt(i)/(R_gas*T1_opt(i));
    viscosity1_opt(i)=viscosity0*((T1_opt(i)/T0)^(3/2))*((T0+S)/(T1_opt(i)+S));

    T2_opt(i)=T01_opt(i)*(1-(gamma-1)/2*Mr_opt(i)^2*(c2_U_opt^2));
    T02_opt(i)=T01_opt(i);
    M2_opt(i)=c2_opt/sqrt(gamma*R_gas*T2_opt(i));
    
    T03_opt(i)=T01_opt(i)*(1-(gamma-1)*Mr_opt(i)^2*psi_opt);
    T3_opt(i)=T03_opt(i)- T01_opt(i)*((gamma-1)/2*Mr_opt(i)^2*(c1_U_opt^2));
    M3_rel_opt(i)=w3_opt/sqrt(gamma*R_gas*T3_opt(i));

    %To calculate the isentropic temperatures and to get the pressures we
    %need to check if the Reynolds correction should be applied. The
    %problem is coupled so we must iterate to get the proper result,
    %initially we will assume that Re=10^5
    max_iter=20; 
    tol=1e-3;
    Re_s_ratio(i)=1; 
    Re_r_ratio(i)=1;

    for j=1:max_iter
        
        loss_coeff_r.Re_corr_opt(i)=loss_coeff_r.Hb_corr_opt*Re_r_ratio(i);
        loss_coeff_s.Re_corr_opt(i)=loss_coeff_s.Hb_corr_opt*Re_s_ratio(i);

        T2s_opt(i)=T2_opt(i)-loss_coeff_s.Re_corr_opt(i)*c2_opt^2/(2*cp);
        p2_opt(i)=p01_opt(i)*(T2s_opt(i)/T01_opt(i))^(gamma/(gamma-1));
        p02_opt(i)=p2_opt(i)/(T2_opt(i)/T02_opt(i))^(gamma/(gamma-1));
        rho2_opt(i)=p2_opt(i)/(R_gas*T2_opt(i)); 
        viscosity2_opt(i)=viscosity0*((T2_opt(i)/T0)^(3/2))*((T0+S)/(T2_opt(i)+S));

        T3s_opt(i)=T3_opt(i)-loss_coeff_r.Re_corr_opt(i)*w3_opt^2/(2*cp);
        p3_opt(i)=p02_opt(i)*(T3s_opt(i)/T02_opt(i))^(gamma/(gamma-1)); 
        p03_opt(i)=p3_opt(i)/(T3_opt(i)/T03_opt(i))^(gamma/(gamma-1));
        rho3_opt(i)=p3_opt(i)/(R_gas*T3_opt(i));
        viscosity3_opt(i)=viscosity0*((T3_opt(i)/T0)^(3/2))*((T0+S)/(T3_opt(i)+S));


        Area_opt(i)=m_dot/(rho2_opt(i)*phi_opt*U_init);
        H_opt(i)=Area_opt(i)/(2*pi*rm_init);
        
        s_over_b_s_opt(i)=0.4/(cos(alpha_2_opt)^2*(tan(alpha_2_opt)+tan(alpha_1_opt)));
        s_over_b_r_opt(i)=0.4/(cos(beta_3_opt)^2*(tan(beta_2_opt)+tan(beta_3_opt)));
        
        b_opt(i)=H_opt(i)/H_over_b;
        s_s_opt(i)=s_over_b_s_opt(i)*b_opt(i);
        s_r_opt(i)=s_over_b_r_opt(i)*b_opt(i);
        
        Dh_s_opt(i)=(2*s_s_opt(i)*H_opt(i)*cos(alpha_2_opt))/(s_s_opt(i)*cos(alpha_2_opt)+H_opt(i)); 
        Dh_r_opt(i)=(2*s_r_opt(i)*H_opt(i)*cos(beta_3_opt))/(s_r_opt(i)*cos(beta_3_opt)+H_opt(i));

        Re_s_opt(i)=rho2_opt(i)*c2_opt*Dh_s_opt(i)/viscosity2_opt(i);
        Re_r_opt(i)=rho3_opt(i)*w3_opt*Dh_r_opt(i)/viscosity3_opt(i);
        Re_s_ratio_new(i)=(1e5/Re_s_opt(i))^0.25;
        Re_r_ratio_new(i)=(1e5/Re_r_opt(i))^0.25;
   
        if abs(Re_s_ratio_new(i)/Re_s_ratio(i)-1)<=tol && abs(Re_r_ratio_new(i)/Re_r_ratio(i)-1)<=tol
            Re_s_ratio(i)= Re_s_ratio_new(i);
            Re_r_ratio(i)= Re_r_ratio_new(i);
            loss_coeff_r.Re_corr_opt(i)=loss_coeff_r.Hb_corr_opt*Re_r_ratio(i);
            loss_coeff_s.Re_corr_opt(i)=loss_coeff_s.Hb_corr_opt*Re_s_ratio(i);
            
            T2s_opt(i)=T2_opt(i)-loss_coeff_s.Re_corr_opt(i)*c2_opt^2/(2*cp);
            p2_opt(i)=p01_opt(i)*(T2s_opt(i)/T01_opt(i))^(gamma/(gamma-1));
            p02_opt(i)=p2_opt(i)/(T2_opt(i)/T02_opt(i))^(gamma/(gamma-1));
            rho2_opt(i)=p2_opt(i)/(R_gas*T2_opt(i)); 
            viscosity2_opt(i)=viscosity0*((T2_opt(i)/T0)^(3/2))*((T0+S)/(T2_opt(i)+S));
    
            T3s_opt(i)=T3_opt(i)-loss_coeff_r.Re_corr_opt(i)*w3_opt^2/(2*cp);
            p3_opt(i)=p02_opt(i)*(T3s_opt(i)/T02_opt(i))^(gamma/(gamma-1)); 
            p03_opt(i)=p3_opt(i)/(T3_opt(i)/T03_opt(i))^(gamma/(gamma-1));
            rho3_opt(i)=p3_opt(i)/(R_gas*T3_opt(i));
            viscosity3_opt(i)=viscosity0*((T3_opt(i)/T0)^(3/2))*((T0+S)/(T3_opt(i)+S));

            break;
        end
        Re_s_ratio(i)= Re_s_ratio_new(i);
        Re_r_ratio(i)= Re_r_ratio_new(i);
    end
    A1_opt(i)=m_dot/(rho1_opt(i)*phi_opt*U_init);
    H1_opt(i)=A1_opt(i)/(2*pi*rm_init);
    b1_opt(i)=H1_opt(i)/H_over_b;
    r1_max(i)=rm_init+H1_opt(i)/2;

    A2_opt(i)=m_dot/(rho2_opt(i)*phi_opt*U_init);
    H2_opt(i)=A2_opt(i)/(2*pi*rm_init);
    b2_opt(i)=H2_opt(i)/H_over_b;
    r2_max(i)=rm_init+H2_opt(i)/2;

    A3_opt(i)=m_dot/(rho3_opt(i)*phi_opt*U_init);
    H3_opt(i)=A3_opt(i)/(2*pi*rm_init);
    b3_opt(i)=H3_opt(i)/H_over_b;
    r3_max(i)=rm_init+H3_opt(i)/2;

    bracket_s_opt(i)=T01_opt(i)/T2_opt(i)*(loss_coeff_s.Re_corr_opt(i)/(1-loss_coeff_s.Re_corr_opt(i)))*c2_U_opt^2;
    bracket_r_opt(i)=(T01_opt(i)/T3_opt(i))*(loss_coeff_r.Re_corr_opt(i)/(1-loss_coeff_r.Re_corr_opt(i)))*w3_U_opt^2;
    denominator_opt(i)=1+(1/(2*psi_opt))*(T03_opt(i)/T01_opt(i))*(bracket_s_opt(i)+bracket_r_opt(i));
    efficiency_opt(i)=1/denominator_opt(i);
    pol_efficiency(i)=(gamma/(gamma-1))*log(T03_opt(i)/T01_opt(i))/log(p03_opt(i)/p01_opt(i));

end

%Calculations for variables needed for lab 2 or asked for this one

p_ratio=p03_opt(end)/p01_opt(1);
T03s=T01_opt(1)*p_ratio^((gamma-1)/gamma);
machine_tt_eff=(T01_opt(1)-T03_opt(end))/(T01_opt(1)-T03s);
machine_pol_eff=(gamma/(gamma-1))*log(T03_opt(end)/T01_opt(1))/log(p_ratio);
machine_overall=machine_pol_eff*machine_tt_eff;

Re_1s=rho1_opt(1)*c1_opt*b1_opt(1)/viscosity1_opt(1);
Re_1r=rho2_opt(1)*c2_opt*b2_opt(1)/viscosity2_opt(1);

M2_rel_1=w2_opt/sqrt(gamma*R_gas*T2_opt(1));

%Lastly we will plot the annulus lines with the leading and trailing edges

H_stages=[H1_opt;H2_opt;H3_opt];

figure;
hold on;

x_annulus=zeros(1, 2*stages+1);
x_annulus(1)=0;
y_hub=zeros(1, 2*stages+1);
y_hub(1)=rm_init-0.5*H1_opt(1);
y_tip=zeros(1, 2*stages + 1);
y_tip(1)=rm_init+0.5*H1_opt(1);
actual_x = 0;

for i = 1:stages
    H_in=H1_opt(i);
    H_mid=H2_opt(i);
    H_out=H3_opt(i);
    b_s = H_in/4;
    b_r = H_mid/4;

    %For the stator
    x_le_s=actual_x;
    x_te_s=x_le_s+b_s;

    rhub_le_s=rm_init-0.5*H_in;   
    rtip_le_s=rm_init+0.5*H_in;
    rhub_te_s=rm_init-0.5*H_mid;  
    rtip_te_s=rm_init+0.5*H_mid;
    position_x_stator=[x_le_s x_te_s x_te_s x_le_s x_le_s];
    position_y_stator=[rhub_le_s rhub_te_s rtip_te_s rtip_le_s rhub_le_s];

    fill(position_x_stator, position_y_stator,[88 88 88]/255);

    gap=0.4*b_s;
    x_le_r=x_te_s+gap;

    % For the rotor
    x_te_r=x_le_r+b_r;

    rhub_le_r=rhub_te_s;  
    rtip_le_r=rtip_te_s;
    rhub_te_r=rm_init-0.5*H_out; 
    rtip_te_r=rm_init+0.5*H_out;
    position_x_rotor=[x_le_r x_te_r x_te_r x_le_r x_le_r];
    position_y_rotor=[rhub_le_r rhub_te_r rtip_te_r rtip_le_r rhub_le_r];

    fill(position_x_rotor, position_y_rotor, "red");

    %For the next stage
    actual_x=x_te_r+0.4*b_r;
    x_annulus(2*i)=(x_te_s+x_le_r)/2;
    y_hub(2*i)=rhub_le_r;
    y_tip(2*i)=rtip_le_r;
    x_annulus(2*i+1)=x_te_r;
    y_hub(2*i+1)=rhub_te_r;
    y_tip(2*i+1)=rtip_te_r;
end

% Líneas de anillo y radio medio
plot(x_annulus, y_hub, "k-", "LineWidth", 3);
plot(x_annulus, y_tip, "k-", 'LineWidth', 3);
plot([min(x_annulus) max(x_annulus)], [rm_init rm_init],"--","LineWidth",1.5);
title("Stators and rotors");
xlabel("Width (m)");
ylabel("Radius (m)");
axis equal;
legend("Stator 1","Rotor 1", "Casing Upper", "Casing Lower", "Mean Radius",  "Location","northeast");


%We will plot the results obtained

fprintf("---------------------- DESIGN ----------------------\n");
fprintf("Flow coefficient=%.4f | Loading coefficient=%.4f | \n" + ...
    "Mean radius=%.3f m | R=%.3f | Mean velocity=%.1f m/s | Stages=%d\n", ...
    phi_opt, psi_opt, rm_init, R, U_init, stages);
fprintf("Angles (deg): alpha_1=%.1f | alpha_2=%.1f | beta_2=%.1f | beta_3=%.1f\n", ...
    alpha_1_opt*180/pi, alpha_2_opt*180/pi, beta_2_opt*180/pi, beta_3_opt*180/pi);
fprintf("Velocities (m/s): c_x=%.1f | c_1=%.1f | c_2=%.1f | w_2=%.1f | w_3=%.1f\n", ...
    phi_opt*U_init, c1_opt, c2_opt, w2_opt, w3_opt);
fprintf("-----------------------------------------------------\n\n");


for i = 1:stages
    fprintf("---------------------- STAGE %d ----------------------\n", i);
    fprintf("Mach number: M1=%.4f | M2=%.4f | M3_rel=%.4f\n", Mr_opt(i), M2_opt(i), M3_rel_opt(i));
    fprintf("Height (m): H1=%.4f | H2=%.4f | H3=%.4f\n", H1_opt(i), H2_opt(i), H3_opt(i));
    fprintf("Max radius (m): r_1=%.4f | r_2=%.4f | r_3=%.4f\n", r1_max(i), r2_max(i), r3_max(i));
    fprintf("Efficiencies: total to total=%.4f | polytropic=%.4f\n", efficiency_opt(i), pol_efficiency(i));
    fprintf("Loss (Re): Stator=%.4e | Rotor=%.4e\n", ...
        loss_coeff_s.Re_corr_opt(i), loss_coeff_r.Re_corr_opt(i));

    fprintf("T (K): T1=%.1f | T2=%.1f | T3=%.1f\n", T1_opt(i), T2_opt(i), T3_opt(i));
    fprintf("p (Pa): p1=%.3e | p2=%.3e | p3=%.3e\n", p1_opt(i), p2_opt(i), p3_opt(i));
    fprintf("-----------------------------------------------------\n\n");
end
fprintf("---------------------- MACHINE EFFICIENCIES ----------------------\n");
fprintf("Overall=%.4f | Total to total=%.4f | Polytropic=%.4f\n",machine_overall,machine_tt_eff, machine_pol_eff);

fprintf("---------------------- FOR LAB 2 ----------------------\n");
fprintf("Stator flow angles (deg): alpha_1=%.1f | alpha_2=%.1f\n", alpha_1_opt*180/pi, alpha_2_opt*180/pi);
fprintf("Rotor relative flow angles (deg): beta_2=%.1f | beta_3=%.1f\n", beta_2_opt*180/pi, beta_3_opt*180/pi);
fprintf("Stator Mach numbers: M_1=%.5f | M_2=%.5f\n", Mr_opt(1), M2_opt(1));
fprintf("Rotor relative Mach numbers: M_1_rel=%.5f | M_2_rel=%.5f\n", M2_rel_1, M3_rel_opt(1));
fprintf("Pith to axial-chord ratios: s_s/b_s=%.3f | s_r/b_r=%.1f\n", s_over_b_s_opt(1), s_over_b_r_opt(1));
fprintf("Reynolds number based on inlet conditions and axial chord: " + ...
    "Re_1s=%.3f | Re_1r=%.1f\n", Re_1s, Re_1r(1));
fprintf("-----------------------------------------------------\n\n");


%Function for velocity triangles. Done as a function as it may be needed
%for future labs

function [alpha_1,alpha_2,beta_2,beta_3,epsilon_s,epsilon_r,c1_U,c1,c2_U,c2,w2_U,w2,w3_U,w3]=turbine_triangle(PHI,PSI,R,U)
    alpha_1=atan((PSI/2+R-1)./PHI);
    alpha_2=atan((PSI/2-R+1)./PHI);
    beta_2=atan((PSI/2-R)./PHI);
    beta_3=atan((PSI/2+R)./PHI);
    epsilon_s=alpha_1+alpha_2;
    epsilon_r=beta_2+beta_3;
    c1_U=hypot(PHI,(PSI/2+R-1));
    c1=c1_U*U;
    c2_U=hypot(PHI,(PSI/2+1-R));
    c2=c2_U*U;
    w2_U=hypot(PHI,(PSI/2-R));
    w2=w2_U*U;
    w3_U=hypot(PHI,(PSI/2+R));
    w3=w3_U*U;
end