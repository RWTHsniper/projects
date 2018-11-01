err1 = importdata('p2_6/error1.txt');
err2 = importdata('p2_6/error2.txt');
err10 = importdata('p2_6/error10.txt');
err100 = importdata('p2_6/error100.txt');
%err1000 = importdata('p2_6/error3.txt');


r_time1 = importdata('p2_6/r_time1.txt');
r_time2 = importdata('p2_6/r_time2.txt');
r_time10 = importdata('p2_6/r_time10.txt');
r_time100 = importdata('p2_6/r_time100.txt');
%r_time1000 = importdata('p2_6/r_time3.txt');

semilogy(r_time1,err1,'*-','markersize',10);
hold on;
semilogy(r_time10,err10,'*-','markersize',10);
semilogy(r_time100,err100,'*-','markersize',10);
%semilogy(r_time1000,err1000,'*-','markersize',10);

legend('\nu1=1','\nu1=10','\nu1=100');
title('Error graph (n=7, \nu2=1)');
xlabel('Runtime (s)');
ylabel('Relative error');
hold off;


