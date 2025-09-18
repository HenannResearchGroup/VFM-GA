function [v] = linear_inter(lb,lv,ub,uv,b)
    v = (b-lb)/(ub-lb)*(uv-lv)+lv;
end

