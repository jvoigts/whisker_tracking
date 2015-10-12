
%%

c=0;




vid = VideoReader(fullfile('example_data/example_video.avi'));
Nframes=vid.NumberOfFrames


mousepresent=zeros(1,Nframes);
mousepresent_graded=mousepresent;
mousepresent_graded_mean=mousepresent;
%% detect mouse, this reduces the N of frames to track later

ll=linspace(0,1,50);


skipN=20; % skip every N frames for detecting presence
Nprint=10;
tic;
c=0;
for i=1:Nframes
    if mod(i,skipN)==1
        c=c+1;
        I=read(vid,i); I=(double(I(1:4:end,1:4:end,1))-0)/20;
        
        if i==1
            lastI=I;
        end;
        
        Id=abs(I-lastI);
        lastI=I;
        
        ii=quantile(Id(:),.995) ;
        mousepresent(i:i+skipN)= ii>0.5;
        mousepresent_graded(i:i+skipN)= ii;
        mousepresent_graded_mean(i:i+skipN)= mean(Id(:));
        
        
        if mod(c,Nprint)==0
            t=toc/Nprint; % time in s per frame
            fprintf('%d/%d frames (%d%%) (%d fps)\n',i,Nframes,floor((i*100)/Nframes),round(1/t));
            clf; hold off;
            subplot(2,1,1);
            imagesc(Id); hold on;
            colormap(gray);
            if mousepresent(i)
                %plot((ll+10).*20,-h+30,'r');
                plot([0 size(Id,1)], [0 size(Id,2)],'r')
            else
                %  plot((ll+10).*20,-h+30,'b');
            end;
            title(num2str(ii));
            subplot(2,1,2);
            if i>5000
                hold on;
                
                plot(mousepresent(end-5000:end),'r');
                plot(mousepresent_graded_mean(end-5000:end));
            end;
            drawnow;
            tic;
        end;
    end;
end;

%crop back
mousepresent=mousepresent(1:Nframes);
mousepresent_graded=mousepresent_graded(1:Nframes);
mousepresent_graded_mean=mousepresent_graded_mean(1:Nframes);

disp('done checking mouse presence');
fprintf('%d frames (%.1f%% of input) to track further\n',sum(mousepresent),100*mean(mousepresent));

%% track animal position frame by frame

tracknose.mpos=NaN(1,Nframes); % mean position

tracknose.lpos=tracknose.mpos; % left border of animal
tracknose.rpos=tracknose.mpos; % right "

tracknose.lnose=tracknose.mpos; % left nose position, like lpos but lower treshold
tracknose.rnose=tracknose.mpos; % right "

tracknose.retractgap=tracknose.mpos; % platform position
tracknose.basegap=tracknose.mpos;    % base platform position


skipn=2; % optionally track only every 2nd frame here
ifplot=1; % swicth whether to plot images or not

trackframes=find(mousepresent_graded>1); % only further track frames with a mouse present
c=0;  se = strel('ball',3,3);
skipi=[0:skipn-1];
for fnum=trackframes(1:skipn:end)
    c=c+1;
    if mod(c,20)==0
        if ~ifplot
            clf;
            subplot(211); hold on;
            conf=mousepresent_graded_mean;
            plot(tracknose.mpos.*(conf>0),'k');
            plot(tracknose.lpos.*(conf>0),'g--');
            plot(tracknose.rpos.*(conf>0),'r--');
            plot(tracknose.lnose.*(conf>0),'g');
            plot(tracknose.rnose.*(conf>0),'r');
            %   plot(conf.*60)
            xlim([-1000 0]+fnum);
            drawnow;
        end;
        fprintf('%d/%d frames (%d%%/ %d%%) \n',fnum,Nframes,round((fnum/Nframes)*100),round(((c*skipn)/numel(trackframes))*100));
    end;
    
    I=read(vid,fnum);
    Icrop=I(170:450,10:420,1);
    
    Icrop_nogray=Icrop;
    
    
    widthscale=[1:size(Icrop_nogray,2)];
    
    Icrop_nogray=imdilate(Icrop_nogray,se); % pre-process image for tracking rough features
    
    ff=[ones(1,5).*1,ones(1,5).*-1 ];
    Ig=max(conv2(double(Icrop_nogray),-ff,'same'),0);
    Ig(Ig>350)=0;
    Ig(1:70,:)=0; Ig(:,1:300)=0; Ig(:,400:end)=0;
    
    [~,gp]= max(mean(Ig));
    tracknose.retractgap(fnum+skipi)=gp;
    
    
    Ig=max(conv2(double(Icrop_nogray),ff,'same'),0);
    Ig(Ig>120)=0; Ig(:,1:5)=0; Ig(:,250:end)=0;
    [~,bp]= max(mean(Ig));
    tracknose.basegap(fnum+skipi)=bp;
    Icrop_nogray(:,gp-2:end)=120;
    
    % now get rid of dark stripe on left (mobile) side in case the base
    % platform is too close
    
    dark_left = max(find(max(Icrop_nogray)<80));
    if numel(dark_left)==0; dark_left=0; end;
    dark_left=min(dark_left,150);
    Icrop_nogray(:,1:dark_left+15)=120;
    Icrop_nogray(1:40,1:dark_left+100)=120;
    
    %  Icrop_nogray(Icrop_nogray>10 & Icrop_nogray<50) = 120;
    
    tracknose.mpos(fnum+skipi)=sum((mean((Icrop_nogray<10))./sum(mean((Icrop_nogray<10)))).*widthscale);
    try
        tracknose.lpos(fnum+skipi)=  min(find(mean((Icrop_nogray<10))>.1));
    end;
    try
        tracknose.rpos(fnum+skipi)=  max(find(mean((Icrop_nogray<10))>.1));
    end;
    
    try
        tracknose.lnose(fnum+skipi)=min(find((mean(imerode(imdilate(Icrop_nogray(70:end,:),se),se)<50))>.01));
    end;
    try
        tracknose.rnose(fnum+skipi)=max(find((mean(imerode(Icrop_nogray(40:end,:),se)<20))>.01));
    end;
    
    %  video_run_cnn_01();
    
    if ifplot==1
        clf;
        subplot(211); hold on;
        conf=mousepresent_graded_mean;
        plot(tracknose.mpos.*(conf>0),'k');
        plot(tracknose.lpos.*(conf>0),'g--');
        plot(tracknose.rpos.*(conf>0),'r--');
        plot(tracknose.lnose.*(conf>0),'g');
        plot(tracknose.rnose.*(conf>0),'r');
        xlim([-1000 0]+fnum);
        
        subplot(212);
        hold off; plot(0); hold on;
        imagesc(Icrop_nogray); daspect([1 1 1]);
        
        % subplot(212);
        % imagesc(iout); daspect([1 1 1]);
        %   imagesc(Icrop_nogray);
        %
        hold on; colormap gray;
        
        plot(tracknose.mpos(fnum).*[1 1],[40 50],'k-');
        text(tracknose.mpos(fnum), 45,'mean pos','color','k');
        
        plot([1 1].*tracknose.lpos(fnum),[0 50],'g-');
        text(tracknose.lpos(fnum), 25,'l pos','color','g');
        
        plot([1 1].*tracknose.lnose(fnum),[100 250],'g-');
        text(tracknose.lpos(fnum), 150,'l nose','color','g');
        
        plot([1 1].*tracknose.rpos(fnum),[0 50],'r-');
        text(tracknose.rpos(fnum), 25,'r pos','color','r');
        
        plot([1 1].*tracknose.rnose(fnum),[100 250],'r-');
        text(tracknose.rnose(fnum), 150,'r nose','color','r');
        
        plot([1 1].*tracknose.basegap(fnum),[240 280],'b--');
        text(tracknose.basegap(fnum), 270,'base gap','color','b');
        
        plot([1 1].*tracknose.retractgap(fnum),[240 280],'b--');
        text(tracknose.retractgap(fnum), 270,'retracting  gap','color','b');
        
        drawnow;
    end;
end;
disp('done tracking');
%% assemble results

% find retracting gap pos
ll=[1:size(Icrop,2)];
h=histc(tracknose.retractgap,ll);
h(400:end)=0;
f=normpdf([-5:5],0,2); f=f./sum(f);
h=conv2(h,f,'same');
h(h<100)=0;
clf; hold on;
plot(h,'k','LineWidth',2);
dh=conv(h,[-1 1],'same');
ddh=conv(h,[-1 2 -1],'same');
plot(dh,'r');
plot(ddh,'b');

p=find( sign(dh)~=sign([0 dh(1:end-1)]) );
p(ddh(p)<1)=[]


if numel(p)==1
    disp('only one gap pos detected, adding 2nd by blind assumption');
    p(2)=p(1)+7;
end;

plot(p,0,'rs');

[epochs.gappos]=p(1)-120;

% and by how much we retract
[epochs.retractiondist]=p(1)-p(2);
%% process rest of events

epochs.nosedist_track=epochs.nosedist.*0;

clf; hold on;
conf=tracknose.mouse_confidence;
conf(isnan(conf))=0;
plot((tracknose.mpos.*(conf>0))-80,'k');
plot((tracknose.lpos.*(conf>0))-80,'b--');
plot((tracknose.rpos.*(conf>0))-80,'r--');
plot((tracknose.lnose.*(conf>0))-80,'b');
plot((tracknose.rnose.*(conf>0)-80),'r');
plot(((conf.*60)-50)-80,'g');

tracknose.maxintensity_clean=tracknose.maxintensity;
tracknose.maxintensity_clean(isnan(tracknose.maxintensity_clean))=0;

plot(tracknose.maxintensity_clean+300,'b');
plot(tracknose.laser_indicator+300,'r');


% clean up target platform tracking
tracknose.retractgap_clean=tracknose.retractgap;
tracknose.retractgap_clean(isnan(tracknose.retractgap))=epochs.gappos; % set to base where not known
tracknose.retractgap_clean(tracknose.retractgap_clean<epochs.gappos)=epochs.gappos;

wsize=6;
tracknose.retractgap_median=tracknose.retractgap_clean;
for i=wsize+1:numel(tracknose.retractgap)-wsize-1
    tracknose.retractgap_median(i) = median(tracknose.retractgap([-wsize wsize]+i));
end;
f=normpdf([-10:10],0,5); f=f./sum(f);
tracknose.retractgap_median=conv( tracknose.retractgap_median,f,'same');

plot(tracknose.retractgap+300,'g');
plot(tracknose.retractgap_median+300,'k');
%plot((tracknose.retractgap.*0)+epochs.gappos+300,'g--');
plot((tracknose.retractgap.*0)+epochs.gappos-epochs.retractiondist+300,'g--');

plot(tracknose.basegap+300,'g');

%plot(epochs.nosedist_display(1,:),'c');
%xlim([95000  100000]+3500);

% replace downward nose zeros with pos

t=find((tracknose.lnose < 1) .* (conf>0)   );
tracknose.lnose_fix=tracknose.lnose;
tracknose.rnose_fix(t)=tracknose.lpos(t)-15;


% fix transient tracking breakdowns
lnose_diff = [0 abs(diff(tracknose.lnose))];
brkdwn=find((lnose_diff > 50) .* (conf>0)   );
if numel(brkdwn) >0
    plot(brkdwn,100,'ks');
    
    for i=1:numel(brkdwn)-1
        t=brkdwn(i)-1:brkdwn(i+1)+1;
        if numel(t) >20
            t=t(1:20);
        end;
        d=tracknose.lnose(t(1)-0)-tracknose.lpos(t(1)-0);
        dnext=tracknose.lnose(t(1)+1)-tracknose.lpos(t(1)+1);
        if abs(dnext)>40
            tracknose.lnose_fix(t)=tracknose.lpos(t)+d;
        end;
    end;
    plot((tracknose.lnose_fix.*(conf>0))-80,'b');
end;
% on the rightward ones well just use rpos when nose is poast the cutoff
% for the gap

t=find((tracknose.rnose-80 > 290) .* (conf>0)   );
tracknose.rnose_fix=tracknose.rnose;
tracknose.rnose_fix(t)=tracknose.rpos(t)+15;
plot((tracknose.rnose_fix.*(conf>0))-80,'k');


%identify absent periods
f=normpdf([-500:500],0,200); f=f./sum(f);
conf_sm=conv(conf,f,'same');
conf_sm_tr=conf_sm>0.1;
plot((conf_sm_tr).*60-50,'k');

onsets=find(diff(conf_sm_tr)==1);
offsets=find(diff(conf_sm_tr)==-1);
onsetpos=[]; offsetpos=[];

mousetimes=find(conf>0.05);

trs=100;
for i=1:numel(onsets)
    onsets(i)= min(mousetimes(mousetimes>onsets(i)));
    onsetpos(i)=tracknose.mpos( onsets(i) );
end;


for i=1:numel(offsets)
    offsets(i) = max(mousetimes(mousetimes<offsets(i)));
    offsetpos(i)=tracknose.mpos( offsets(i) );
end;

offsets=[1 offsets];
offsetpos=[onsetpos(1) offsetpos];

onsets=[onsets numel(conf)];
onsetpos=[onsetpos offsetpos(end)];

plot(offsets,offsetpos,'bs');
plot(onsets,onsetpos,'rs');

trs=200;
np_lo=-100;
np_hi=300;

epochs.gapdist=epochs.nosedist.*NaN;

for i=1:numel(onsets)-1
    t=[offsets(i):onsets(i)];
    tt=[onsets(i):offsets(i+1)];
    if  offsetpos(i) < trs % towzrds retracting platform
        epochs.nosedist_track(1,t)=np_lo;
        epochs.nosedist_track(1,tt)=tracknose.rnose_fix(tt)-80;
    else    % return
        epochs.nosedist_track(1,t)=np_hi;
        epochs.nosedist_track(1,tt)=tracknose.lnose_fix(tt)-80;
    end;
    
    epochs.nosedist_track(2,:)=200;
    
    
    % per crossing, determine base gap dist
    %[~,epochs.gapdist]=0;
    
    d=tracknose.basegap(tt);
    ll=[1:max(d)];
    h=histc(d,ll);
    [~,gp]=max(h);
    epochs.gapdist([tt(1)-100:tt(end)+100])=gp-120;
    
    
    % and where we retracted the target platform
    if  offsetpos(i) < trs % towzrds retracting platform
        retracted=tracknose.retractgap_median(tt)-120>mean([epochs.gappos epochs.gappos-epochs.retractiondist]);
        retracted(1)=0;
        retracted(end)=0;
        
        on=find(diff(retracted)==1);
        % offsets=find(diff(retracted)==-1);
        
        for j=1:numel(on)
            epochs.data_frames(end+1,1)=on(j)+tt(1);
            epochs.data_frames(end,2)=on(j)+5+tt(1);
            epochs.data_frames(end,3)=8;
            
            plot([0 15]+on(j)+tt(1),[1 1].*700,'r','LineWidth',15);
        end;
        
        
    end;
    
    
end;
% find laser times
laser_on=(tracknose.maxintensity_clean>200) | (tracknose.laser_indicator>30);
laser_on=conv(double(laser_on),ones(1,200),'same')>0;



laser_on(1)=0;laser_on(end)=0;
on=find((diff(laser_on)==1) );
off=find(diff(laser_on)==-1);
for j=1:numel(on)
    epochs.data_frames(end+1,1)=on(j);
    epochs.data_frames(end,2)=off(j);
    epochs.data_frames(end,3)=12;
    
    plot([on(j) off(j)],[1 1].*100,'b','LineWidth',15); drawnow;
end;


%xlim([0 2000]+327800);

plot(epochs.nosedist_track(1,:),'k','LineWidth',2);
plot(epochs.nosedist(1,:),'r','LineWidth',2);

u=find(~isnan(epochs.gapdist));
epochs.gapdist_display=interp1(u,epochs.gapdist(u), [1:numel(epochs.gapdist)]);

epochs.nosedist=epochs.nosedist_track-40;
epochs.nosedist_display=epochs.nosedist_track-40;



