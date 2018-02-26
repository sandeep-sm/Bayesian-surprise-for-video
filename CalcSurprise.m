m=128;
minx = 0;
maxx = 10;
N = 100001;
x = linspace(minx,maxx,N);
surprise = zeros(m,m,nframes);
resimg = zeros(m,m,nframes);
decay_factor = 0.5;
for i = 1:1:nframes
    resimg(:,:,i) = imresize(outimg(:,:,i), [m m]);
end

for i = 1:1:m
    for j = 1:1:m
        data = [];
        for k = 1:1:nframes
            data = [data,resimg(i,j,k)];
        end
        if resimg(i,j,1) ~= 0
            smod = newsm(decay_factor,resimg(i,j,1),resimg(i,j,1),x);
        else
            smod = newsm(decay_factor,0.001,0.001,x);
        end
        
        smod.options.debug = 0;
        smod.options.factordecay = 'yes';
        smod.options.graph = 'no';
        % run the model on the data you have. 
        smod = runsm(data,smod,x);
        surprise(i,j,:) = smod.surprise; 
    end
end


