% example video tracking setup
%
% Used here for gap-crossing experiments with an AVT Pike F-032B 
% using the image aq. toolbox
% http://www.mathworks.com/hardware-support/allied-vision-technologies.html
%
% 2014 Jakob Voigts (jvoigts@mit.edu)

imaqreset

vid = videoinput('avtmatlabadaptor64_r2009b', 1, 'F7M0_Mono8_640x480');
src = getselectedsource(vid);

vid.FramesPerTrigger = Inf;

vid.ROIPosition = [125 0 441 469];

src.Timebase = 0; % select microseconds
src.Shutter = 170; % in \mus

diskLogger = VideoWriter('E:\filename.avi', 'Uncompressed AVI');
% E is a raid 0 - you'll need either an SSD or a good raid 0 to handle the
% data rate comong off the camera in this configuration.

vid.LoggingMode = 'disk';
vid.DiskLogger = diskLogger;

imaqmem( 600*1024^3); % crashes after 16gig of data without this (why?)
%% optional, run preview
preview(vid);

%%
stoppreview(vid);

%%
start(vid)
%%
stop(vid)
