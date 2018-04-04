#/bin/bash!
FNAME=uc.txt

rm $FNAME

# uncertainties to run over
# UC=( -1 -2 -3 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 -19 0 1 2 3 9 10 11 12 13 14 15 16 17 18 19 20 )
# UC=( -3 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18  3 9 10 11 12 13 14 15 16 17 18 )
# UC=( -1 -2 -19 -21 0 1 2 19 21 20 )
UC=( 22, 23 )
echo Running over Uncertainties
for i in "${UC[@]}"; do
    echo $i
    printf "%s\n" $i >> $FNAME
done

# first run over the data to produce raw histos
# then projections

#./analysis 0 -1 0 0 1

# make spectra, dPhi (default) histograms
./analysis 1 -1 0 0 1
./analysis 3 -1 0 0 1

# next. run over the MC's with whatever uncertainties (include nominal)

for ucf in "${UC[@]}"; do
    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf MC RAW
    ./analysis 0 -1 0 0 0 0 $ucf

    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf MC SPECT, CF
    ./analysis 1 -1 0 0 0 0 $ucf

    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf MC DPHI, CF
    ./analysis 3 -1 0 0 0 0 $ucf
    
    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf DATA SPECT UNFOLDING
    ./analysis 2 -1 0 0 1 0 $ucf
    
    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf DATA DPHI UNFOLDING
    ./analysis 4 -1 0 0 1 0 $ucf

done

rm $FNAME
