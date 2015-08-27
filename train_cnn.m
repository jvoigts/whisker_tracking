
% Z:\ephys_ntsr1\convnet_nosetracking

%% make training set
%load('C:\Users\jvoigts\Dropbox\wh_tracking\convnet_nosetracking\trainingset2.mat');


train_dir='Z:\ephys_ntsr1\convnet_whtracking\training_jvs';

tset=[];
tset.in=[];
tset.out=[];
c=0;

locations_pos=[];
locations_neg=[];
for i=[1 4 6]
    c=c+1;
    I  = imread(fullfile(train_dir,[num2str(i),'a.png']));
    imagesc(I); drawnow;
    tset.in(:,:,c)=imread(fullfile(train_dir,[num2str(i),'.png']));
    tset.out(:,:,c)=I>5;
    
    locations_pos=[locations_pos; find(I(:))];
    
    r=zeros(size(I));
    r(randi(numel(I),numel(find(I(:)))*30,1))=1;
    r(1:20,:)=0; r(end-20:end,:)=0; r(:,1:20)=0; r(:,end-20:end)=0;
    r(find(I(:)))=0;
    locations_neg=[locations_neg; find(r(:))];
end;


locations=[locations_pos; locations_neg];
locations_labels=[locations_pos>0; locations_neg.*.0];

nlocations=numel(locations);
kTrainNum = nlocations


inradius=5;

clear TrainX TrainY TrainY_im;
c=0;
for inflatetraining=1:3
for x=1:nlocations
    
    
    c=c+1;
    ii=x;%ceil(nlocations*rand);
    [i,j,k]=ind2sub(size(I),locations(ii));
    
    tile=double(( tset.in(i-inradius:i+inradius,...
        j-inradius:j+inradius,k)./255)-0.5);
    
    for i=1:randi(4)-1
        tile=rot90(tile);
    end;
    
    TrainX(:,:,x)=tile+randn*.5;
    
    
    if ~locations_labels(x)
        TrainY(x,:)=[0 1]';
    else
        TrainY(x,:)=[1 0]';
    end;
    c=c+1;
end;
end;

%% run


%funtype = 'gpu';
%funtype = 'cpu';
funtype = 'matlab';

disp(funtype);


kSampleDim = ndims(TrainX);
kXSize = size(TrainX);
kXSize(kSampleDim) = [];
if (kSampleDim == 3)
    kXSize(3) = 1;
end;
kWorkspaceFolder = './workspace';
if (~exist(kWorkspaceFolder, 'dir'))
    mkdir(kWorkspaceFolder);
end;


kOutputs = size(TrainY, 2);
train_x = single(TrainX(:, :, 1:kTrainNum));
train_y = single(TrainY(1:kTrainNum, :));

kTestNum = nlocations;
test_x = single(TrainX(:, :, 1:kTestNum));
test_y = single(TrainY(1:kTestNum, :));

clear params;
params.epochs = 1;
params.alpha = 0.2;
% this is the parameter for invariant backpropagation
% keep it 0 for standard backpropagation
params.beta = 0;
params.momentum = 0.9;
params.lossfun = 'logreg';
params.shuffle = 1;
params.seed = 0;
dropout = 0;

% norm_x = squeeze(mean(sqrt(sum(sum(train_x.^2))), kSampleDim));

% !!! IMPORTANT NOTICES FOR GPU VERSION !!!
% Outputmaps number should be divisible on 16
% For speed use only the default value of batchsize = 128

% This structure gives pretty good results on MNIST after just several epochs

layers = {
    struct('type', 'i', 'mapsize', kXSize(1:2), 'outputmaps', kXSize(3))
    % remove the following layer in the Matlab version - it is not implemented there
    %  struct('type', 'j', 'mapsize', [28 28], 'shift', [1 1], ...
    %         'scale', [1.40 1.40], 'angle', 0.10, 'defval', 0)
    % struct('type', 'c', 'filtersize', [5 5], 'outputmaps', 8)
    %struct('type', 's', 'scale', [2 2], 'function', 'max', 'stride', [2 2])
    %struct('type', 'c', 'filtersize', [5 5], 'outputmaps', 64, 'padding', [2 2])
    %struct('type', 's', 'scale', [3 3], 'function', 'max', 'stride', [2 2])
    struct('type', 'f', 'length', 8, 'dropout', dropout)
    struct('type', 'f', 'length', 4, 'dropout', dropout)
    struct('type', 'f', 'length', kOutputs, 'function', 'soft')
    };

%weights = single(genweights(layers, params, funtype));
EpochNum = 400;
errors = [];%zeros(EpochNum, 1);

 clf;
 
for i = 1 : EpochNum
    disp(['Epoch: ' num2str((i-1) * params.epochs + 1)])
    [weights, trainerr] = cnntrain(layers, weights, params, train_x, train_y, funtype);
    disp([num2str(mean(trainerr(:, 1))) ' loss']);
    [err, bad, pred] = cnntest(layers, weights, params, test_x, test_y, funtype);
    disp([num2str(err*100) '% error']);
    errors(i) = err;
   % params.alpha = params.alpha * 0.999;
   % params.beta = params.beta * 0.999;
   %
    subplot(2,2,1);
    hold off;
    plot(errors); drawnow;
    
    if mod(i,15)==0
        % ----
        imstack=[];
        for k=1%:size(tset.in,3)
            uim = tset.in(:,:,k);
            % make tiles to feed into CNN
            x=0;
            scalefactor=1;
            clear pos imstack;
            isteps=inradius+10:scalefactor:size(uim,1)-inradius-10;
            jsteps = inradius+10:scalefactor:size(uim,2)-inradius-10;
            for i=isteps
                for j=jsteps
                    x=x+1;
                    imstack(:,:,x)=double((uim(i-inradius:i+inradius,j-inradius:j+inradius)./255)-0.5);
                    
                end;
            end;
            
            pred = cnnclassify(layers, weights, params, imstack, funtype);
            iout=reshape(pred(:,2),numel(jsteps),numel(isteps));
            
            subplot(2,2,[3 4]);
            hold off;
            imagesc(((iout))); colormap(gray); daspect([1 1 1]);
            drawnow;
        end;
        % ----
    end;
end;
disp('Done!');

%% save conv net
% save('cnn_8_3_15.mat','layers','weights', 'params','inradius');

%%
 load('cnn_8_3_15.mat')
 
funtype = 'matlab';
%% run on image

imstack=[];
for k=3%:size(tset.in,3)
    uim = tset.in(:,:,k);
    % make tiles to feed into CNN
    x=0;
    scalefactor=1;
    clear pos imstack;
    isteps=inradius+10:scalefactor:size(uim,1)-inradius-10;
    jsteps = inradius+10:scalefactor:size(uim,2)-inradius-10;
    for i=isteps
        for j=jsteps
            x=x+1;
            imstack(:,:,x)=double((uim(i-inradius:i+inradius,j-inradius:j+inradius)./255)-0.5);
            
        end;
    end;
    
    pred = cnnclassify(layers, weights, params, imstack, funtype);
    iout=reshape(pred(:,2),numel(jsteps),numel(isteps));
    
    clf; hold on;
    imagesc(iout); colormap(gray); daspect([1 1 1]);
    drawnow;
end;
%[err, bad, pred] = cnntest(layers, weights, params, test_x, test_y, funtype);
