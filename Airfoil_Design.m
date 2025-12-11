clc
clear
close all

t_over_c=0.2;
a=0.4;
true_chord=1;
t=t_over_c*true_chord;
pitch_to_chord_s=0.966;
pitch_to_chord_r=1.3;

%Load the thickness excel
data=readmatrix("Profile_Thickness_Distributions.xlsx","Sheet","A3K7");
x_thickness=data(:,1)/100;
y_thickness=data(:,2)/100;
valid=~isnan(x_thickness)&~isnan(y_thickness);
x_thickness=x_thickness(valid);
y_thickness=y_thickness(valid);


[x_thickness,iu]=unique(max(0,min(1,x_thickness)),"stable");
y_thickness=y_thickness(iu);
t_fun=@(xx) interp1(x_thickness,y_thickness,xx,"pchip","extrap");


%For the stator
alpha1=0.385882669398074*180/pi;
alpha2=-1.176927716058535*180/pi;
incidence=0;
metal1_s=alpha1;
flow_deflect_s=alpha1-alpha2; 
deviation_s=1.5/100*flow_deflect_s;
metal2_s=deviation_s+alpha2;

%Añadir stager en alpha1 y restarlo en alpha2

stagger_s=0.5*((metal2_s)+(metal1_s));

metal1_s_new=metal1_s+stagger_s;
metal2_s_new=metal2_s-stagger_s;

A_s=[0,0,1,0,0,0;
    0,1,0,0,0,0;
    0,0,0,1,1,1;
    0,0,0,2,1,0;
    a^2,a,1,-(a^2),-a,-1;
    2*a,1,0,-2*a,-1,0];
b_s=[0;
    tand(metal1_s-stagger_s);
    0;
    tand(metal2_s-stagger_s);
    0;
    0];
x_s=A_s\b_s;

a1_s=x_s(1); 
b1_s=x_s(2); 
c1_s=x_s(3);
a2_s=x_s(4); 
b2_s=x_s(5); 
c2_s=x_s(6);

x_values_1_s=linspace(0,a,200);
y_values_1_s=a1_s*x_values_1_s.^2+b1_s*x_values_1_s+c1_s;
x_values_2_s=linspace(a,1,200);
y_values_2_s=a2_s*x_values_2_s.^2 + b2_s*x_values_2_s + c2_s;
x_s=[x_values_1_s(1:end-1),x_values_2_s];
y_s=[y_values_1_s(1:end-1), y_values_2_s];

figure; 
hold on; 
grid on;
plot([x_values_1_s x_values_2_s],[y_values_1_s y_values_2_s],"LineWidth",1.8);
title("Stator camber line");
xlim([0 1]);
ylim([0 1]);

%To add the thickness
x_clip=max(0,min(1,x_s));
t_s=t_fun(x_clip); 

dy_s=gradient(y_s, x_s);
phi_s=atan(dy_s);  

xu_s=x_s-0.5.*t_s.*sin(phi_s);
yu_s=y_s+0.5.*t_s.*cos(phi_s);

xl_s=x_s+0.5.*t_s.*sin(phi_s);
yl_s=y_s-0.5.*t_s.*cos(phi_s);

% To plot the rotated airfoil:
Rs=[cosd(stagger_s) -sind(stagger_s);
    sind(stagger_s)  cosd(stagger_s)];

US=Rs*[xu_s;yu_s];  
LS=Rs*[xl_s;yl_s];   

xu_s_rot=US(1,:); 
yu_s_rot=US(2,:);

xl_s_rot=LS(1,:); 
yl_s_rot=LS(2,:);

figure; 
hold on; 
axis equal; 
grid on;
plot(xu_s_rot,yu_s_rot, "b", "LineWidth", 1.3);
plot(xl_s_rot, yl_s_rot, "r", "LineWidth", 1.3);
title("Stator profile rotated");
xlabel("x/c"); 
ylabel("y/c");

name="Stator.g6s";
mises_file(name,xu_s_rot,yu_s_rot,xl_s_rot,yl_s_rot,alpha1,alpha2,pitch_to_chord_s)

%For the rotor
beta2=-0.093476781158589*180/pi;
beta3=-1.239399769330577*180/pi;
metal2_r=beta2;
flow_deflect_r=beta2-beta3; 
deviation_r=1.5/100*flow_deflect_r;
metal3_r=deviation_r+beta3;

stagger_r=0.5*((metal3_r)+(metal2_r));

A_r=[0,0,1,0,0,0;
    0,1,0,0,0,0;
    0,0,0,1,1,1;
    0,0,0,2,1,0;
    a^2,a,1,-(a^2),-a,-1;
    2*a,1,0,-2*a,-1,0];
b_r=[0;
    tand(metal2_r-stagger_r);
    0;
    tand(metal3_r-stagger_r);
    0;
    0];
x_r=A_r\b_r;

a1_r=x_r(1); 
b1_r=x_r(2); 
c1_r=x_r(3);
a2_r=x_r(4); 
b2_r=x_r(5); 
c2_r=x_r(6);

x_values_1_r=linspace(0,a,200);
y_values_1_r=a1_r*x_values_1_r.^2+b1_r*x_values_1_r+c1_r;
x_values_2_r=linspace(a,1,200);
y_values_2_r=a2_r*x_values_2_r.^2+b2_r*x_values_2_r+c2_r;

figure; 

hold on; 
grid on;
plot([x_values_1_r x_values_2_r],[y_values_1_r y_values_2_r],"LineWidth",1.8);
title("Rotor camber line");
xlim([0 1]);
ylim([0 1]);

x_r = [x_values_1_r(1:end-1), x_values_2_r];
y_r = [y_values_1_r(1:end-1), y_values_2_r];
x_clip_r=max(0,min(1,x_r));
t_r=t_fun(x_clip_r);

dy_r=gradient(y_r, x_r);
phi_r=atan(dy_r);

xu_r=x_r-0.5.*t_r.*sin(phi_r);
yu_r=y_r+0.5.* t_r.*cos(phi_r);

xl_r=x_r+0.5.*t_r.*sin(phi_r);
yl_r=y_r-0.5.*t_r.*cos(phi_r);

% To plot the rotated airfoil:
Rr=[cosd(stagger_r) -sind(stagger_r);
    sind(stagger_r)  cosd(stagger_r)];

UR=Rr*[xu_r;yu_r];  
LR=Rr*[xl_r;yl_r];   

xu_r_rot=UR(1,:); 
yu_r_rot=UR(2,:);

xl_r_rot=LR(1,:); 
yl_r_rot=LR(2,:);

figure; 
hold on; 
axis equal; 
grid on;
plot(xu_r_rot,yu_r_rot, "b", "LineWidth", 1.3);
plot(xl_r_rot, yl_r_rot, "r", "LineWidth", 1.3);
title("Rotor profile rotated");
xlabel("x/c"); 
ylabel("y/c");

name="Rotor.g6r";
mises_file(name,xu_r_rot,yu_r_rot,xl_r_rot,yl_r_rot,beta2,beta3,pitch_to_chord_r)


function mises_file(name,xu,yu,xl,yl,alpha_in,alpha_out,pitch)
    SINL=tand(alpha_in);
    SOUT=tand(alpha_out);
    CHINL=0.7;
    CHOUT=0.7;
    PITCH=pitch;
    xu=xu(end:-1:1);
    yu=yu(end:-1:1);
    xl(1)=[];
    yl(1)=[];
    X=[xu,xl];
    Y=[yu,yl];

    fid=fopen(name, "w");
    fprintf(fid,"%f BLADE", name);
    fprintf(fid,"%.6f %.6f %.6f %.6f %.6f\n", SINL, SOUT, CHINL, CHOUT, PITCH);
    for i = 1:length(X)
        fprintf(fid, "%.6f %.6f\n", X(i), Y(i));
    end
    fclose(fid);
end



