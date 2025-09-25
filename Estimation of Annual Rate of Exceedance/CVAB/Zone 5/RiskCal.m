
%==========================================================================================
% Fragility estimation, hazard analysis, Risk calculation, and R_avl calculation from IDA
% Ketan Gupta, IIT Roorkee, 2024
% ketan_g@eq.iitr.ac.in
% Geometric mean of intensities; Risk integral truncated for lower and upper bounds
%==========================================================================================

clear; clc;

%%%%%%%%%% INPUTS %%%%%%%%%%

MCE = 0.6;
median_Op = 0.28; std_Op = 0.39;
median_CP = 0.57; std_CP = 0.29;

% Integral lower bound -
low_int_IO = 0.24;
low_int_CP = 0.24;

% Integral upper bound -
high_int_IO = 0.90;
high_int_CP = 0.90;

hazard = [0.24 0.0137; 0.3 0.00444; 0.4 0.00211; 0.48 0.00103; 0.6 0.00040; 0.75 0.0002; 0.9 0.0001];

%%%%%%%%%% END OF INPUTS %%%%%%%%%%

phi_MCE_Op = 1/std_Op*log(MCE/median_Op);
phi_MCE_CP = 1/std_CP*log(MCE/median_CP);

P_MCE_Op = normcdf(phi_MCE_Op)*100;     % P[Op|MCE] for Operational level
P_MCE_CP = normcdf(phi_MCE_CP)*100;     % P[CP|MCE] for Collapse level
fprintf ("Probability of Exceeding Operational given MCE is %6.3f %%\n",P_MCE_Op);
fprintf ("Probability of Exceeding Collapse given MCE is %6.3f %%\n",P_MCE_CP);

% Discretization -
disStep = 100;
dis_IO = (high_int_IO-low_int_IO)/disStep;
dis_CP = (high_int_CP-low_int_CP)/disStep;

% Hazard Curve -
spec_acc_fr_IO_Risk = low_int_IO:dis_IO:high_int_IO;                                    % Hazard lower and upper limits
spec_acc_fr_CP_Risk = low_int_CP:dis_CP:high_int_CP;

for i = 1:size(spec_acc_fr_CP_Risk,2)
    % hazrad_mod(i) = interp1(hazard(:,1),hazard(:,2),spec_acc_fr_CP_Risk(i));
    hazrad_mod(i) = exp(interp1(hazard(:,1),log(hazard(:,2)),spec_acc_fr_CP_Risk(i)));
end


figure
loglog(spec_acc_fr_CP_Risk,hazrad_mod,'k','LineWidth',0.75)
xlabel('Spectral Acceleration, Sa_g_e_o (g)')
ylabel('Annual Rate of Exceedance')
res = 1200;
box on
ax = gca;
ax.LineWidth = 1.2;
set(gca, 'FontName', 'Times New Roman')
x0=500;
y0=300;
width=220;
height=180;
set(gcf,'units','points','position',[x0,y0,width,height]);
print('Hazard Curve.png','-dpng',['-r' num2str(res)]);

% Fragility function -
fragility_Op = makedist('Lognormal','mu',median_Op,'sigma',std_Op);
fragility_CP = makedist('Lognormal','mu',median_CP,'sigma',std_CP);

% Fragility Curve -
frag_high = norminv(0.9999,median_CP, std_CP);                                         % Upper limit for fragility plot
spec_acc_fr = 0:dis_CP:frag_high;
fragility_Op_plot = normcdf((log(spec_acc_fr/median_Op))/std_Op);
fragility_CP_plot = normcdf((log(spec_acc_fr/median_CP))/std_CP);
figure
plot(spec_acc_fr,fragility_Op_plot, '-r', 'linewidth', 1)
hold on
plot(spec_acc_fr,fragility_CP_plot, '-b', 'linewidth', 1)
hold on
grid on
%title('Fragility Curve')
xlabel('Spectral Acceleration, Sa(T1) (g)')
ylabel('Probability of Exceedance')
xlim([0 frag_high])
yticks(0:0.25:1);
%legend('Fragility Curve');
med_IO = round(median_Op,2);
med_CP = round(median_CP,2);
std_IO = round(std_Op,2);
std_CP = round(std_CP,2);
dim = [0.55, 0.2, 0.3, 0.15];
str = {['\mu_I_O = ',num2str(med_IO),', ','\beta_I_O = ',num2str(std_IO)],['\mu_C_P = ',num2str(med_CP),', ','\beta_C_P = ',num2str(std_CP)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize', 6, 'FontName', 'Times New Roman', 'BackgroundColor', 'w');
res = 1200;
box on
ax = gca;
ax.LineWidth = 1.2;
set(gca, 'FontName', 'Times New Roman')
x0=500;
y0=300;
width=220;
height=180;
set(gcf,'units','points','position',[x0,y0,width,height]);

% Modification in fragility for area = 1
fragility_IO_mod = truncate(fragility_Op,low_int_IO,high_int_IO);
fragility_CP_mod = truncate(fragility_CP,low_int_CP,high_int_CP);


plot(spec_acc_fr,cdf(fragility_IO_mod,spec_acc_fr),':r','LineWidth',1.5)
plot(spec_acc_fr,cdf(fragility_CP_mod,spec_acc_fr),':b','LineWidth',1.5)
print('Fragility Curve.png','-dpng',['-r' num2str(res)]);

% Fragility for Risk integration -
ii = 0;
for i = low_int_IO:dis_IO:high_int_IO
    ii = ii + 1;
    fragility_Risk_IO(ii) = cdf(fragility_IO_mod,i);
end
ii = 0;
for i = low_int_CP:dis_CP:high_int_CP
    ii = ii + 1;
    fragility_Risk_CP(ii) = cdf(fragility_CP_mod,i);
end

% Risk -
m_IO = size(spec_acc_fr_IO_Risk,2);
for i = 2:m_IO-1
    del_s = spec_acc_fr_IO_Risk(i) - spec_acc_fr_IO_Risk(i-1);
    del_F = fragility_Risk_IO(i) - fragility_Risk_IO(i-1);
    mi = log(hazrad_mod(i)/hazrad_mod(i-1))/del_s;
    ai = hazrad_mod(i)*(1-exp(mi * del_s));
    bi = hazrad_mod(i)/del_s*(exp(mi*del_s)*(del_s-1/mi)+1/mi);
    lambda_IO(i) = fragility_Risk_IO(i)*ai - del_F*bi;
end

sum3_IO = 0;
size3 = size(lambda_IO,2);
for i = 1:size3
    sum3_IO = sum3_IO + lambda_IO(i);
end
mean_ann_R_IO = sum3_IO;
ann_prob_IO = 1-exp(-1*mean_ann_R_IO);
Mean_ann_R_IO = mean_ann_R_IO*100;
Ann_prob_IO = ann_prob_IO*100;
fprintf ("Mean Annual Rate of Exceeding IO is %6.3f %%\n",Mean_ann_R_IO);

m_CP = size(spec_acc_fr_CP_Risk,2);
for i = 2:m_CP-1
    del_s = spec_acc_fr_CP_Risk(i) - spec_acc_fr_CP_Risk(i-1);
    del_F = fragility_Risk_CP(i) - fragility_Risk_CP(i-1);
    mi = log(hazrad_mod(i)/hazrad_mod(i-1))/del_s;
    ai = hazrad_mod(i)*(1-exp(mi * del_s));
    bi = hazrad_mod(i)/del_s*(exp(mi*del_s)*(del_s-1/mi)+1/mi);
    lambda_CP(i) = fragility_Risk_CP(i)*ai - del_F*bi;
end

sum3_CP = 0;
size3 = size(lambda_CP,2);
for i = 1:size3
    sum3_CP = sum3_CP + lambda_CP(i);
end
mean_ann_R_CP = sum3_CP;
ann_prob_CP = 1-exp(-1*mean_ann_R_CP);
Mean_ann_R_CP = mean_ann_R_CP*100;
Ann_prob_CP = ann_prob_CP*100;
fprintf ("Mean Annual Rate of Exceeding CP is %6.3f %%\n",Mean_ann_R_CP);

% Save results to text file -
fid = fopen('Fragility & Risk Results.txt', 'w');
fprintf(fid,"Median for IO = %5.3f g ; Std. Dev. for IO = %5.3f; \nMedian for CP = %5.3f g ; " + ...
    "Std. Dev. for CP = %5.3f; \n\nP[IO|MCE] = %5.2f %%\nP[CP|MCE] = %5.2f %%\n\nMean Annual Rate of Exceeding IO   = %6.3f %%" + ...
    "\nMean Annual Rate of Exceeding CP   = %6.3f %%", median_Op, std_Op, median_CP, std_CP, P_MCE_Op, P_MCE_CP, Mean_ann_R_IO, Mean_ann_R_CP);
fclose(fid);