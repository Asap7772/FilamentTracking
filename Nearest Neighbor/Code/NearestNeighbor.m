%----------------------------------------------------------------------------------------------
%A method to ensure that each frame has the same number of points
%Author: Anikait Singh
%Date: July 19th, 2017
%Version Number: 1
%----------------------------------------------------------------------------------------------
function [ pos ] = NearestNeighbor( I_skel, frames)
%% Finding Frame with Maximum Number of Points
tic
[row, ~] = find(I_skel(:,:,1));
[maxPoints, ~] = size(row);
index = 1;
for i = 2:frames
    [row, ~] = find(I_skel(:,:,i));
    [numPoints,~] = size(row);
    if(maxPoints > numPoints)
        maxPoints = numPoints;
        index = i;
    end
end
toc
%% Link similar points together
% Create the new matrix
tic
disp(frames)
disp(maxPoints)
disp(index)
pos = zeros(maxPoints, 2, frames);
pos(pos == 0) = NaN;
toc
%% From max frame backwards to first frame
tic
for i = index-1:-1:1
    
    %case when max frame is the first frame
    [rowPrev, colPrev] = find(I_skel(:,:,i+1));
    [rowCurr, colCurr] = find(I_skel(:,:,i));
    
    for x = 1: size(rowCurr)
        close = 1;
        X = [rowCurr(x), colCurr(x);rowPrev(1),colPrev(1)]; %matrix of two points
        d = pdist(X,'euclidean'); % finds distance
        for y = 2: size(rowPrev)
            X = [rowCurr(x), colCurr(x);rowPrev(y),colPrev(y)];
            currDist = pdist(X,'euclidean');
            if(currDist < d)
                d = currDist;
                close = y;
            end
        end
        
        %TODO what if another point is there?
        pos(close,:,i) = [colCurr(x), rowCurr(x)];
    end
end

%% From max frame to end
for i = index+1:frames
    
    %case when max frame is the first frame
    [rowPrev, colPrev] = find(I_skel(:,:,i-1));
    [rowCurr, colCurr] = find(I_skel(:,:,i));
    
    for x = 1: size(rowCurr)
        close = 1;
        X = [rowCurr(x), colCurr(x);rowPrev(1),colPrev(1)]; %matrix of two points
        d = pdist(X,'euclidean'); % finds distance
        for y = 2: size(rowPrev)
            X = [rowCurr(x), colCurr(x);rowPrev(y),colPrev(y)];
            currDist = pdist(X,'euclidean');
            if(currDist < d)
                d = currDist;
                close = y;
            end
        end
        
        %TODO what if another point is there?
        pos(close,:,i) = [colCurr(x), rowCurr(x)];
    end
end
toc
%% Ensure that max frame's points show up
[row, col] = find(I_skel(:,:,index));
pos(1:size(row),1,index) = col(1:size(row));
pos(1:size(col),2,index) = row(1:size(col));
%% Labeling points as NAN if they disapear
pos(pos ==0) = nan;
%% Plotting
for i = 1:frames
    if(mod(i,10) == 0)
        subplot(2,5,10)
        plot(pos(:,1,i), pos(:,2,i), '*')
        axis([0 350 0 100]);
        title(strcat("Frame", num2str(i)));
        if(not(i == 100))
            figure
        end
    else
        subplot(2,5,mod(i,10))
        plot(pos(:,1,i), pos(:,2,i), '*')
        axis([0 350 0 100]);
        title(strcat("Frame", num2str(i)));
    end
end