%% 2.2 f) Hough Peaks for a 3D Hough Matrix
function peaks = houghpeaks3d(H,numpeaks)
    peaks = zeros(numpeaks,3);
    for i = 1:numpeaks
        [~,idx] = max(H(:));
        [idx1,idx2,idx3] = ind2sub(size(H),idx);
        peaks(i,:) = [idx1,idx2,idx3];
        H(idx1-40:idx1+40,idx2-40:idx2+40,:) = 0;
    end
end