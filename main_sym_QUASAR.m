close all
clear

%Variable and simulation parameter
udm_factor_scale = 6000;
r=1:5;
f=[90, 50, 20]./udm_factor_scale;
s=[1 3 5];
nTI = 13;
snr=[100,20,10,5];
step   = .3;
t   = [.04:step:(step*(nTI-1)+.04)];
nr=100;
T1b=1.6;
alpha = 2*pi/360*35;
look_locker_scale = exp(t.*(-log(cos(alpha))/(t(2)-t(1))));
delay = [0,0.3,0.6,0.9];
p = 0.1;

%Preallocation of results

aif                  = zeros(length(f),length(r),length(s),length(snr),length(delay),length(t));
r_true               = zeros(length(f),length(r),length(s),length(snr),length(delay),length(t));
rm_true              = zeros(length(f),length(r),length(s),length(snr),length(delay),length(t));
DM_true              = zeros(length(f),length(r),length(s),length(snr),length(delay),length(t));
sigma                = zeros(length(f),length(r),length(s),length(snr),length(delay),length(t));
DM_noise             = zeros(length(f),length(r),length(s),length(snr),length(delay),length(t),nr);



%% - Simulation iteration


for fi=1:length(f)                          % Perfusion
    for ri=1:5                                % Residue function
        for si=1:length(s)                  % AIF shape
            for snri=1:length(snr)          % SNR 
                for di=1:length(delay)      % Delay 
                    %MC - Simulation
                    [aif(fi,ri,si,snri,di,:), ...
                        r_true(fi,ri,si,snri,di,:), ...
                        rm_true(fi,ri,si,snri,di,:), ...
                        DM_true(fi,ri,si,snri,di,:), ...
                        DM_noise(fi,ri,si,snri,di,:,:), ...
                        sigma(fi,ri,si,snri,di,:)]=symASL_QUASAR(t,f(fi),ri,s(si),delay(di),snr(snri),nr);
                end
            end
        end
    end
end


%% Results
for fi=3 % 1 Hyper - 2 GM - 3 WM
    figure(fi)
    j=1;
    for ri=[1 4 5 3] % plot Residue functions
        subplot(2,2,j)
        for snri=4:-1:2 % plot SNR
            plot_ci(t,squeeze(DM_noise(fi,ri,1,snri,1,:,:)),5,50,95,[0 snri/4 0],[0 0 0])
            hold on
        end
        xlim([t(1) t(end)])
        ylim([-.001 .003])
        j=j+1;
    end
end

