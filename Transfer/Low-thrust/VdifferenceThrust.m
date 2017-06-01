ReverseTime = xlsread('C:\Users\Bart\Desktop\DSE\DSE_Moon_TBD\Transfer\Data\SpiralBackTimeDatastep001.xlsx');
SpiralOut = xlsread('C:\Users\Bart\Desktop\DSE\DSE_Moon_TBD\Transfer\Data\SpiralOutDatastep001.xlsx');

Thrust = ReverseTime(:,1);
VReverseTime = sqrt(ReverseTime(:,6).^2 + ReverseTime(:,7).^2);
VSpiralOut = sqrt(SpiralOut(:,6).^2 + SpiralOut(:,7).^2);


Vdifference = VReverseTime-VSpiralOut;
time = ReverseTime(:,2)+SpiralOut(:,2);
time = time/(60*60*24);
plot(Thrust,abs(Vdifference))
xlabel('Thrust [N]')
ylabel('V_{difference} [m/s]')
title('Difference in velocity at L1 between phase 1 and phase 2')

