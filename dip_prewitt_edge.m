%% 1.2.1) Prewitt Edge Detector function
function edge_img = dip_prewitt_edge(img,thresh)
    Gpx = [-1 0 1; -1 0 1; -1 0 1]/6;
    Gpy = [1 1 1; 0 0 0 ;-1 -1 -1]/6;
    Gx = conv2(img,Gpx,"same");
    Gy = conv2(img,Gpy,"same");
    edge_img = sqrt(Gx.^2+Gy.^2);
    edge_img(edge_img<thresh) = 0; edge_img(edge_img>=thresh) = 1;
end