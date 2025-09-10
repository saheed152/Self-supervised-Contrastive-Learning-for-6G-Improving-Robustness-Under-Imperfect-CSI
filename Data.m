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
addpath 'C:\Users\sahee\Desktop\Matlab Files\benchmarks\AltMinAlg\Narrowband\PE-AltMin'
Ns = 2; %number of data streams to be transmitted
NRF =2;  %This is Nrf for one user. actual NRF of the system = NRF*user
User=5;
SNR_dB = -15:5:30;
SNR = 10.^(SNR_dB./10);    
smax = length(SNR);% enable the parallel
In_PE=0;
%WRF=1
power = .1;
Nr=16;
Nt=144;
%load("awgn_noise.mat")
%load('H_new_36_144_5_200.mat')
%load('all_u5_ns2.mat')
%load("ICC_DL_10000.mat");
load("ICC_channel_Nr16_1000real.mat");
realization_brute=1;
noise_pow=0.1;
r_sel = 4;
all_Comb_Fopt = zeros(realization_brute,Nt,User*r_sel);
all_Comb_sv = zeros(realization_brute,User*r_sel,User*r_sel);
selected_index_real = zeros(realization_brute,User*Ns);
%sum_max = zeros(realization_brute);
%time = zeros(realization_brute);
for z=1:realization_brute
    z
    H=H_new(:,:,:,z);
    Comb_Fopt = zeros(Nt,User*r_sel);
    Comb_Wopt = zeros(Nr,User*r_sel);
    Fopt_corr = zeros(Nt,User*r_sel);
    Comb_sv = zeros(User*r_sel);
    tStart_brute=0;
    tStart_brute=tic;
    for i=1:User
        [U,S,V] = svd(H(:,:,i));
        Comb_Fopt(:,(i-1)*r_sel+1:i*r_sel) = V([1:Nt],[1:r_sel]);
        Comb_Wopt(:,(i-1)*r_sel+1:i*r_sel) = U([1:Nr],[1:r_sel]);
        temp_F = V*S';
        Fopt_corr(:,(i-1)*r_sel+1:i*r_sel) = temp_F([1:Nt],[1:r_sel]); % Getting V multiplied by sigmas
        temp_sv = diag(S);
        Comb_sv((i-1)*r_sel+1:i*r_sel) = temp_sv(1:r_sel);
    end
    all_Comb_Fopt(z,:,:) = Comb_Fopt;
    all_Comb_sv(z,:,:) = Comb_sv;
    u = User;
    r = r_sel;

    A = nchoosek([1:r],Ns);

    a = 1:size(A,1);  % These should contain inded of all combinations of one user's vectors ( rCNs)
    % Put all vectors into cell array
    allVecs = {a,a,a,a,a}; % We have to put as many vectors as there are users.
    sub = cell(1,numel(allVecs));
    [sub{:}] = ndgrid(allVecs{:});
    sub = cellfun(@(x)x(:),sub,'UniformOutput', false);
    allPerms = cell2mat(sub); % These are the indexes of each

    ind = zeros(size(allPerms,1),u*Ns);
    for i = 1:size(allPerms,1)
        ind(i,:) = reshape(A(allPerms(i,:),:)',1,u*Ns);
    end

    % The following loops sets indices for the combined matrix
    for i = 1:u-1
        ind(:,Ns*i+1:Ns*(i+1)) = ind(:,Ns*i+1:Ns*(i+1)) + i*r;
    end

    R_sumo = zeros(length(SNR),1);
    %Max_R = zeros(length(SNR),1); %uncomment for a snr range
    Max_R=0;%for signle SNR
    Min_corr = 10^9;
    max_sv = sum(Comb_sv(ind(1,:)));

    for j = 1:size(ind,1)
        j;
        Fopt = Comb_Fopt(:,ind(j,:));
        Fopt = reshape(Fopt,[Nt,Ns,User]);
        Wopt = Comb_Wopt(:,ind(j,:));
        Wopt = reshape(Wopt,[Nr,Ns,User]);
        %for s = 1:smax
        s=8;%consider snr=8dB
            R_o=0;
            for U=1:User
                In_opt=0;
                Fopt_sum=0;
                for p=1:User
                    if p~=U
                        Fopt_sum=Fopt_sum+Fopt(:,:,p)*Fopt(:,:,p)';
                    end
                end

                Ck=(SNR(s)/User)*Wopt(:,:,U)'*H(:,:,U)*Fopt_sum*H(:,:,U)'*Wopt(:,:,U)+noise_pow.*Wopt(:,:,U)'*Wopt(:,:,U);
                R_o=R_o+log2(det(eye(Nr)+(SNR(s)/User)*Wopt(:,:,U)*inv(Ck)*Wopt(:,:,U)'*H(:,:,U)*Fopt(:,:,U)*Fopt(:,:,U)'*H(:,:,U)'));
            end
            R_sumo=R_o; %Suboptimal
            if ((abs(R_sumo)) > (abs(Max_R)))
                Max_R = abs(R_sumo);
                selected_index=ind(j,:);
            end
        %end

    end
    selected_index_real(z,:)=selected_index;
    sum_max(z)=Max_R;
    tEnd_brute=toc(tStart_brute);
    time(z)=tEnd_brute;
end

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
