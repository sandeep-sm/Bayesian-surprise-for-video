vid = VideoReader('Neovision2-Training-Heli-013.mpg');
n=512;
k = 1;
nframes=40;
img = zeros(n,n,3,nframes);
outimg = zeros(n,n,nframes);
vid.CurrentTime=0;

p = makeGBVSParams;
p.channels = 'CIOM';
p.levels = 3;
p.salmapmaxsize = 40;
p.blurfrac = 0.04;

motinfo = [];

while hasFrame(vid)
    img(:,:,:,k) = imresize(readFrame(vid),[n n]);
    [out{k}, motinfo] = gbvs(img(:,:,:,k) , p , motinfo );
    outimg(:,:,k) = out{k}.master_map_resized;
    if k == nframes
        break;
    end
    k = k+1;
    vid.CurrentTime=vid.CurrentTime+0.15;
end