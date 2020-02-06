m=128;
iterator = 0;
minx = 0;
maxx = 200;
N = 4000001;
x = linspace(minx,maxx,N);
surprise = zeros(m,m,nframes);
resimg = zeros(m,m,nframes);
decay_factor = 1;
for i = 1:1:nframes
    resimg(:,:,i) = imresize(outimg(:,:,i), [m m]);
end

for i = 1:1:m
    for j = 1:1:m
        data = [];
        for k = 1:1:nframes
            data = [data,round(resimg(i,j,k),4)];
        end
        if data(1) ~= 0
            b = data(1);
            smod = newsm(decay_factor,b,b,x);
        else
            b = 0.005;
            smod = newsm(decay_factor,b,b,x);
        end
        
        smod.options.debug = 0;
        smod.options.factordecay = 'no';
        smod.options.graph = 'no';
        % run the model on the data you have. 
        smod = runsm(data,smod,x,iterator);
        iterator = smod.iterator;
        surprise(i,j,:) = smod.surprise; 
    end
end


