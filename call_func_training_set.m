%%%%%%%%%%-----saheed Ullah----%%%%%%
% Modify: 1/24/2024 :- modify the code to generate the data for ICC DL
% paper; Modify it to export the array index for each realizations

%%%%%%% Rafid %%%%%%
%Modify: 1/24/2024: Declared variables to be saved and added command to
%save them

clear,clc
close all
addpath(pwd);
addpath(genpath(pwd));
addpath 'C:\Users\mdsah\Desktop\Xtra Research\research\HybridPrecodingOpt-master\benchmarks\AltMinAlg\Narrowband\PE-AltMin'

load("awgn_noise.mat")
%load('H_new_36_144_5_200.mat')
%load('all_u5_ns2.mat')
load("ICC_channel_Nr16_1000real.mat");
%load("")
realization_brute=10;
noise_pow=0.1;
r_sel = 4;





plot(SNR_dB,sum(sum_max,2)/realization_brute,'b-->','Marker','>','LineWidth',1);
xlabel('SNR [dB]')
ylabel('Spectral efficiency (bits/s/Hz)')
legend('PE-Altmin','MO-AltMin','Sub-optimal','BF_search','Corr_F','Location','NW')

save('allsaved_ICC.mat', 'all_Comb_Fopt', 'all_Comb_sv', 'selected_index_real','time','sum_max');
%
% figure(1);
% fig=figure(1);
% set(fig,'Units','inches');
% screenposition = get(fig,'Position');
% set(fig,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print(fig,'Augmented_Signal','-dpdf','-r0')
