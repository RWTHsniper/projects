
clear;

err1 = importdata('p2_7/error1.txt');
err2 = importdata('p2_7/error2.txt');
err3 = importdata('p2_7/error3.txt');
err7 = importdata('p2_7/error7.txt');
err10 = importdata('p2_7/error10.txt');
err20 = importdata('p2_7/error20.txt');

r_time1 = importdata('p2_7/r_time1.txt');
r_time2 = importdata('p2_7/r_time2.txt');
r_time3 = importdata('p2_7/r_time3.txt');
r_time7 = importdata('p2_7/r_time7.txt');
r_time10 = importdata('p2_7/r_time10.txt');
r_time20 = importdata('p2_7/r_time20.txt');

semilogy(r_time1, err1,'*-','markersize',10);
hold on;
semilogy(r_time2,err2,'d-','markersize',10);
semilogy(r_time3,err3,'d-','markersize',10);
semilogy(r_time7,err7,'d-','markersize',10);
semilogy(r_time10,err10,'d-','markersize',10);
semilogy(r_time20,err20,'d-','markersize',10);

legend('\nu1=1','\nu1=2','\nu1=3','\nu1=7','\nu1=10','\nu1=20');
title('Error graph (n=7,\nu2=1)');
xlabel('Runtime');
ylabel('Relative error');
grid;
hold off;

