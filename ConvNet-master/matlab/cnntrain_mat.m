function [weights, trainerr] = cnntrain_mat(layers, weights, params, train_x, train_y)

params = setparams(params);
layers = cnnsetup(layers, 0);
assert(length(size(train_y)) == 2, 'The label array must have 2 dimensions'); 
train_num = size(train_y, 1);
classes_num = size(train_y, 2);
assert(classes_num == layers{end}.length, 'Labels and last layer must have equal number of classes');
if (params.balance == 1)
  layers{end}.coef = (ones(1, classes_num) ./ mean(train_y, 1)) / classes_num;
elseif (params.balance == 0)
  layers{end}.coef = ones(1, classes_num);  
end;
layers = setweights(layers, weights);

assert(size(train_x, 1) == layers{1}.mapsize(1) && ...
       size(train_x, 2) == layers{1}.mapsize(2), ...
       'Data and the first layer must have equal sizes');
assert(size(train_x, 3) == layers{1}.outputmaps, ...
       'The number of the input data maps must be as specified');
assert(size(train_x, 4) == train_num, ...
       'Data and labels must have equal number of objects');

train_x = normalize(layers{1}, train_x);     
layers{1} = initnorm(layers{1}, train_x);
     
rng(params.seed);
numbatches = ceil(train_num/params.batchsize);
trainerr = zeros(params.numepochs, 2);
for epoch = 1 : params.numepochs  
  if (length(params.beta) == 1)
    beta = params.beta;
  else
    beta = params.beta(epoch);
  end;  
  if (params.shuffle == 0)
    kk = 1:train_num;
  else
    kk = randperm(train_num);
  end;
  for batch = 1 : numbatches
    
    batch_x = train_x(:, :, :, kk((batch-1)*params.batchsize + 1 : min(batch*params.batchsize, train_num)));    
    batch_y = train_y(kk((batch-1)*params.batchsize + 1 : min(batch*params.batchsize, train_num)), :);
    
    % first pass
    layers = initact(layers, batch_x);
    layers = updateweights(layers, params, epoch, 0); % preliminary update
    [layers, pred] = forward(layers, 1);
    %disp(['pred: ' num2str(pred(1,1:5))]);
    
    % second pass
    [layers, loss] = initder(layers, params, batch_y);
    trainerr(epoch, 1) = trainerr(epoch, 1) + loss;    
    %disp(['loss: ' num2str(loss)]);
    layers = backward(layers, params);
    layers = calcweights(layers, 2);
    
    % third pass
    [layers, loss2] = initder2(layers);
    trainerr(epoch, 2) = trainerr(epoch, 2) + loss2;
    %disp(['loss2: ' num2str(loss2)]);      
    if (beta > 0)      
      [layers, pred2] = forward(layers, 3);          
      %disp(['pred2: ' num2str(pred2(1,1:5))]);
      layers = calcweights(layers, 3);      
    end;
    
    layers = updateweights(layers, params, epoch, 1); % final update
     
    if (params.verbose == 2)
      disp(['Epoch: ' num2str(epoch) ', batch: ', num2str(batch)]);
    end;
  end
  trainerr = trainerr / numbatches;
  if (params.verbose == 1)
    disp(['Epoch: ' num2str(epoch)]);
  end;
end
    
weights = getweights(layers); 

end
