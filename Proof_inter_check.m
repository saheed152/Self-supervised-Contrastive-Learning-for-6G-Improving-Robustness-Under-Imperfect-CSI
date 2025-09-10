clear,clc
close all
addpath(pwd);
addpath(genpath(pwd));
addpath 'C:\Users\mdsah\Desktop\Xtra Research\research\HybridPrecodingOpt-master\benchmarks\AltMinAlg\Narrowband\PE-AltMin'
Ns = 2; %number of data streams to be transmitted
NRF =2;  %This is Nrf for one user. actual NRF of the system = NRF*user
User=5;
SNR_dB = -15:5:25;
SNR = 10.^(SNR_dB./10);
smax = length(SNR);% enable the parallel
realiza=1;
In_PE=0;
%WRF=1
power = .1;
Nr=36;
Nt=144;
% n = wgn(Nr,1,power); %AWGN Noise
% n1 = wgn(Ns,1,power);

%load("awgn_noise.mat")


%load('../../datasets/all_u5_ns2_1.mat'); %Load the channnel matrix

load('all_u5_ns2.mat')
%%%PE-AltMin

tStart = tic;
r_sel = 4;
Comb_Fopt = zeros(Nt,User*r_sel);
Comb_Wopt = zeros(Nr,User*r_sel);
Fopt_corr = zeros(Nt,User*r_sel);
Comb_sv = zeros(User*r_sel);

for i=1:User
    [U,S,V] = svd(H(:,:,i));
    Comb_Fopt(:,(i-1)*r_sel+1:i*r_sel) = V([1:Nt],[1:r_sel]);
    Comb_Wopt(:,(i-1)*r_sel+1:i*r_sel) = U([1:Nr],[1:r_sel]);
    temp_F = V*S';
    Fopt_corr(:,(i-1)*r_sel+1:i*r_sel) = temp_F([1:Nt],[1:r_sel]); % Getting V multiplied by sigmas
    Fopt_corr2(:,(i-1)*r_sel+1:i*r_sel) = V([1:Nt],[1:r_sel]);
    temp_sv = diag(S);
    Comb_sv((i-1)*r_sel+1:i*r_sel) = temp_sv(1:r_sel);
end

big_Fopt=Fopt_corr'*Fopt_corr2;
u = User;
r = r_sel;

A = nchoosek([1:r],Ns);

a = 1:size(A,1);  % These should contain inded of all combinations of one user's vectors ( rCNs)
% Put all vectors into cell array
allVecs = {a,a,a,a,a}; % We have to put as many vectors as there are users.
sub = cell(1,numel(allVecs));
[sub{:}] = ndgrid(allVecs{:});
sub = cellfun(@(x)x(:),sub,'UniformOutput', false);
% allPerms is [m x n] matrix of m permutations of n vectors
% m should equal prod(cellfun(@numel,allVecs))
% n should equal numel(allVecs)
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
Max_R = zeros(length(SNR),1);
Min_corr = 10^9;
max_sv = sum(Comb_sv(ind(1,:)));
size_big=r_sel*User;
kount=0;
full_cor_mat = zeros(size_big,size_big)
for j = 1:size(ind,1)
    j
    if (sum(Comb_sv(ind(j,:)))>0.75*max_sv)

        %Fsel_corr = Fopt_corr(:,ind(j,:));
        %temp_diag = Fsel_corr'*Fsel_corr - diag(diag(Fsel_corr'*Fsel_corr));
        temp_corr_sel = sel_matrix(ind(j,:),size_big).*big_Fopt;
        temp_diag = temp_corr_sel - diag(diag(temp_corr_sel));
        Corr_F = norm(temp_diag,'fro');

        if ((Corr_F < Min_corr))
            Min_corr = Corr_F;
            ind_corr_min = j;
            full_cor_mat = temp_diag;
        end
    end

end


Fopt = Comb_Fopt(:,ind(ind_corr_min,:));

% Corr_F = norm(abs(Fopt'*Fopt),'fro')^2 - trace((Fopt'*Fopt).^2);

Fopt = reshape(Fopt,[Nt,Ns,User]);


Wopt = Comb_Wopt(:,ind(ind_corr_min,:));
Wopt = reshape(Wopt,[Nr,Ns,User]);

In_2_corr = zeros(User,1)

for U=1:User
    In_opt=0;
    for p=1:User
        if p~=U
            In_opt=In_opt+norm(Wopt(:,:,U)' * H(:,:,U) * Fopt(:,:,p),'fro')^2;
        end
    end
    In_2_corr(U)=In_opt;

end
tEnd = toc(tStart)

sum(In_2_corr)
Min_corr^2
