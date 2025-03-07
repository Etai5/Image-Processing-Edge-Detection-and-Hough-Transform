%% 2.1 c) Hough Matrix for lines function
function [HM, radius_vec, theta_vec] = dip_hough_lines(BW,radius,theta)
    [M, N] = size(BW);
    radius_vec = -sqrt(M^2+N^2) : radius : sqrt(M^2+N^2);
    theta_vec = -90 : theta : 90;
    HM = zeros(length(radius_vec),length(theta_vec));
    [y_vec, x_vec] = find(BW);
    for i = 1:length(y_vec)
        y = y_vec(i);
        x = x_vec(i);
        for t_index = 1:length(theta_vec)
            r = x*cosd(theta_vec(t_index)) + y*sind(theta_vec(t_index));
            [~,r_index] = min(abs(radius_vec-r));
            HM(r_index,t_index) = HM(r_index,t_index)+1;
        end
    end
end