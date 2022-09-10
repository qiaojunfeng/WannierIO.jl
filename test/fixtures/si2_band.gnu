set style data dots
set nokey
set xrange [0: 5.90043]
set yrange [ -6.79852 :  7.20233]
set arrow from  1.15701,  -6.79852 to  1.15701,   7.20233 nohead
set arrow from  1.56608,  -6.79852 to  1.56608,   7.20233 nohead
set arrow from  2.27460,  -6.79852 to  2.27460,   7.20233 nohead
set arrow from  3.50179,  -6.79852 to  3.50179,   7.20233 nohead
set arrow from  4.50380,  -6.79852 to  4.50380,   7.20233 nohead
set arrow from  5.32193,  -6.79852 to  5.32193,   7.20233 nohead
set xtics ("G"  0.00000,"X"  1.15701,"U"  1.56608,"K"  2.27460,"G"  3.50179,"L"  4.50380,"W"  5.32193,"X"  5.90043)
 plot "si2_band.dat"
