settling_time = 11;
t = time_cpet;
tt = t - settling_time;

MRO2p_0 = bestP(9);
MRO2p_1 = bestP(8);
MRO2p_2 = bestP(7);
MRO2p_3 = bestP(6);
MRO2p_4 = bestP(5);
MRO2p_5 = bestP(4);
MRO2p_6 = bestP(3);
MRO2p_7 = bestP(2); 
MRO2p_8 = bestP(1); 

MRO2 = MRO2p_8*tt.^8 + MRO2p_7*tt.^7 + MRO2p_6*tt.^6 + MRO2p_5*tt.^5 + ...
       MRO2p_4*tt.^4 + MRO2p_3*tt.^3 + MRO2p_2*tt.^2 + ...
       MRO2p_1*tt + MRO2p_0;
MRO2 = MRO2 /60000;
plot(tt, MRO2)