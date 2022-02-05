function [ f1,f2,f3,f4,f5,f6,f7,f8] = getfvd4( ir1,ir2,ir3,ir4,ir5,ir6,ir7,ir8 )

mse=cumsum(ir1.^2)+cumsum(ir2.^2)+cumsum(ir3.^2)+cumsum(ir4.^2)+cumsum(ir5.^2)+...
    cumsum(ir6.^2)+cumsum(ir7.^2)+cumsum(ir8.^2);
f1=cumsum(ir1.^2)./mse;
f2=cumsum(ir2.^2)./mse;
f3=cumsum(ir3.^2)./mse;
f4=cumsum(ir4.^2)./mse;
f5=cumsum(ir5.^2)./mse;

f6=cumsum(ir6.^2)./mse;
f7=cumsum(ir7.^2)./mse;
f8=cumsum(ir8.^2)./mse;




end

