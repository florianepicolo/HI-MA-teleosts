
datapathway="DATA/"
listfile="list_test.txt"
outfilename="a"
argRprog=$outfilename"_chi2.csv"

chmod +x step1.py
chmod +x step2.R

./step1.py -f $datapathway -l $listfile -s $outfilename
R CMD BATCH "--args $argRprog" step2.R
rm step2.Rout