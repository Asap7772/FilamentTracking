%----------------------------------------------------------------------------------------------
%A method to ensure that each frame has the same number of points
%Author: Anikait Singhe
%Date: July 19th, 2017
%Version Number: 1
%----------------------------------------------------------------------------------------------
function [ pos ] = NearestNeighbor_Advanced( I_skel, frames)
%% Finding Frame with Maximum Number of Points
tic
[row, ~] = find(I_skel(:,:,1));
[maxPoints, ~] = size(row);
index = 1;
for ii = 2:frames
    [row, ~] = find(I_skel(:,:,ii));
    [numPoints,~] = size(row);
    if(maxPoints < numPoints)
        maxPoints = numPoints;
        index = ii;
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
for ii = index-1:-1:1
    
    %case when max frame is the first frame
    [rowPrev, colPrev] = find(I_skel(:,:,ii+1));
    [rowCurr, colCurr] = find(I_skel(:,:,ii));
    
    [N, ~] = size(rowCurr);
    for x = 1: N
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
        pos(close,:,ii) = [colCurr(x), rowCurr(x)];
        
        %         disp(strcat(num2str(x), strcat(' index: ', num2str(ii))))
        
        if(not(x >= N) && pos(x,1,ii) == 0)
            pos(x,:,ii) = pos(x,:,ii+1);
        end
    end
end

%% From max frame to end
for ii = index+1:frames
    
    %case when max frame is the first frame
    [rowPrev, colPrev] = find(I_skel(:,:,ii-1));
    [rowCurr, colCurr] = find(I_skel(:,:,ii));
    
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
        pos(close,:,ii) = [colCurr(x), rowCurr(x)];
        if(pos(x,1,ii) == 0)
            pos(x,:,ii) = pos(x,:,ii+1);
        end
    end
end
toc
%% Ensure that max frame's points show up
[row, col] = find(I_skel(:,:,index));
pos(1:size(row),1,index) = col(1:size(row));
pos(1:size(col),2,index) = row(1:size(col));
%% Interpolating points
% From max frame backwards to first frame

tic
for i = index-1:-1:1
    for x = 1: size(rowCurr)
        if(isNAN(pos(x,1,i)))
            pos(x,:,i) = pos(x,:,i+1);
        end
    end
end

% From max frame to end
for i = index+1:frames
    for x = 1: size(rowCurr)
        if(isNAN(pos(x,1,i)))
            pos(x,:,i) = pos(x,:,i+1);
        end
    end
end
toc
%% Plotting
for ii = 1:frames
    if(mod(ii,10) == 0)
        subplot(2,5,10)
        plot(pos(:,1,ii), pos(:,2,ii), '*')
        axis([0 350 0 100]);
        title(strcat("Frame", num2str(ii)));
        if(not(ii == 100))
            figure
        end
    else
        subplot(2,5,mod(ii,10))
        plot(pos(:,1,ii), pos(:,2,ii), '*')
        axis([0 350 0 100]);
        title(strcat("Frame", num2str(ii)));
    end
end