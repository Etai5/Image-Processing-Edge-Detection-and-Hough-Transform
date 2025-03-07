%% 2.2 c) Hough Matrix for circles function
function [HM, radius_vec] = dip_hough_circles(BW,radius,theta)
    [M, N] = size(BW);
    radius_vec = 80 : radius : 100;
    theta_vec = 0 : theta : 360;
    HM = zeros(M,N,length(radius_vec));
    [y_vec, x_vec] = find(BW);
    for i = 1:length(y_vec)
        y = y_vec(i);
        x = x_vec(i);
        for r_index = 1:length(radius_vec)
            for t_index = 1:length(theta_vec)
                a = round(x - radius_vec(r_index)*cosd(theta_vec(t_index)));
                b = round(y - radius_vec(r_index)*sind(theta_vec(t_index)));
                if (a>0 && b>0 && a<=M && b<=N)
                    HM(a,b,r_index) = HM(a,b,r_index)+1;
                end
            end
        end
    end
end