function [Y] = fastWinv(X,scaling,dims)
%Compute W^(-1)*X

%LP-cone
for i=1:dims.l
    Y(i,:) = X(i,:)/scaling.l.wl(i);
end
%SOC-cone

for i=1:length(dims.q) %do this for every cone
    Y_temp = 0;
    X_temp = X(dims.l+sum(dims.q(1:(i-1)))+1:dims.l+sum(dims.q(1:(i))),:);
    for p=1:size(X_temp,2)
        z0 = X_temp(1,p);
        z1 = X_temp(2:end,p);
        zeta = scaling.q(i).q'*z1;
        y = (1/scaling.q(i).eta)*[scaling.q(i).a*z0-zeta;z1+(-z0+zeta/(1+scaling.q(i).a))*scaling.q(i).q]; 
        if p==1
            Y_temp = y;
        else
            Y_temp = [Y_temp y];
        end
    end
    Y = [Y; Y_temp];
end

