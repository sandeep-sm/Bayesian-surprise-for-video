function X = yuckypush(b,X,dim)

if dim == 2
    cols = [];
    for i=1:size(X,2) - 1
        cols = [cols i];
    end
    X = [b'  X(:,cols)]
elseif dim == 1
    rows = [];
    for i=1:size(X,1) -1
        rows = [rows 1];
    end
    X = [b ; X(rows,:)];
else
    error('Diminsion in yuckypush must either be a 1 or 2');
end