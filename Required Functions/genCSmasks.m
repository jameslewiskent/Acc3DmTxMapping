function [masks] = genCSmasks(settings,num_samples_per_acc)
% Idea of this script is that we generate a whole bunch of masks way more
% than we require, then keep the ones closest to the acceleration factor we
% are after, then remove non-unique masks, then cut down to the number of
% requested masks. This is horribly inefficient but it does the job.

sx = settings.sx;
sy = settings.sy;
fovx = settings.fovx;
fovy = settings.fovy;
calib = settings.calib;
ellipse = settings.ellipse;
pp = settings.pp;
accelerations = settings.accelerations;
Mask_N = settings.Mask_N;

dS = 2; % within dS samples of desired accelerations
overshoot_Factor = 1.5; % Go to overshoot_Factor*Mask_N as we will remove duplicates later
tempsettings.sx = sx; tempsettings.sy = sy; tempsettings.fovx = fovx; tempsettings.fovy = fovy; tempsettings.calib = calib; tempsettings.ellipse = ellipse; tempsettings.pp = pp; tempsettings.accelerations = accelerations; tempsettings.Mask_N = Mask_N;

Folder = ['Data',filesep,'Undersampling Masks',filesep,'x',num2str(sx),'y',num2str(sy),'calib',num2str(calib)];
if ~exist(Folder, 'dir')
    mkdir(Folder)
end

if exist(fullfile(Folder,'masks.mat'),'file') == 2
    load(fullfile(Folder,'masks.mat'),'masks','settings')
end

if ~exist('settings','var') || ~isequaln(settings,tempsettings)
    
    if accelerations ~= 0
        num_samples_per_acc = round((sx*sy)./accelerations);
    else
        accelerations = (sx*sy)./num_samples_per_acc; % hence accelerations desired
    end
    masks = ones(sx,sy,size(accelerations,2),ceil(overshoot_Factor*Mask_N)); count = zeros(1,size(accelerations,2));
    for accel_n = 1:size(accelerations,2)
        if accelerations(accel_n) ~= 1
            while any(count(1,accel_n) <= ceil(overshoot_Factor*Mask_N))
                for accel = accelerations(accel_n)*0.75:accelerations(accel_n)*0.01:accelerations(accel_n)*2 % just give wide range
                    accelx = sqrt(accel);
                    accely = accelx;
                    mask = vdPoisMex(sx,sy,fovx,fovy,accelx,accely,calib,ellipse, pp);
                    [minval,minind] = min(abs(sum(mask,1:2) - num_samples_per_acc));
                    if minval < dS % within dS samples of desired accelerations
                        masks(:,:,minind,count(1,minind)+1) = mask;
                        count(1,minind) = count(1,minind) +1;
                    end
                end
            end
        end
        disp(['Acceleration of ',num2str(accelerations(accel_n)),' finished.'])
    end
    masks = masks(:,:,:,1:ceil(overshoot_Factor*Mask_N)); % crop all to same size
    figure('color','w'); imagesc(imtile(masks(:,:,:,1))); title('Example Mask of each Acceleration Factor Provided')
    
    % Check for and remove duplicate masks, then crop down to Mask_N
    Unique_masks = ones([size(masks,1:3),Mask_N]);
    for accel_n = 1:size(accelerations,2)
        clear temp
        if accelerations(accel_n) ~= 1
            temp = reshape(unique(reshape(squeeze(masks(:,:,accel_n,:)),size(masks,1)*size(masks,2),size(masks,4))','rows')',size(masks,1),size(masks,2),1,[]);
            disp(['Number of non-unique masks removed: ', num2str(ceil(overshoot_Factor*Mask_N)-size(temp,4)),'. Masks left: ',num2str(size(temp,4))]);
            Unique_masks(:,:,accel_n,1:Mask_N) = temp(:,:,1,1:Mask_N); % Hopefully we have more unique masks than neccesary...
        end
    end
    masks = Unique_masks;
    
    settings = tempsettings;
    save(fullfile(Folder,'masks.mat'),'masks','settings')
else
    disp('Masks using these settings already exist. Not re-generating.')
end
end