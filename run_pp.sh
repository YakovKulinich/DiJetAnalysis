#/bin/bash!
FNAME=uc.txt

rm $FNAME

# uncertainties to run over
# UC=( -1 -2 -3 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 -19 0 1 2 3 9 10 11 12 13 14 15 16 17 18 19 20 )
UC=( -1 -2 -19 -21 0 1 2 19 20 21 )
echo Running over Uncertainties
for i in "${UC[@]}"; do
    echo $i
    printf "%s\n" $i >> $FNAME
done

# first run over the data to produce raw histos
# then projections

#./analysis 0 -1 0 0 1
#./analysis 1 -1 0 0 1

# next. run over the MC's with whatever uncertainties (include nominal)

for ucf in "${UC[@]}"; do
    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf MC RAW
    # ./analysis 0 -1 0 0 0 0 $ucf

    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf MC SPECT
    ./analysis 1 -1 0 0 0 0 $ucf

    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf MC SPECT UNFOLDING
    ./analysis 2 -1 0 0 0 0 $ucf

    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf MC DPHI
    ./analysis 3 -1 0 0 0 0 $ucf
    
    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf DATA SPECT UNFOLDING
    ./analysis 2 -1 0 0 1 0 $ucf
    
    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf DATA DPHI
    ./analysis 3 -1 0 0 1 0 $ucf

    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf DATA DPHI UNFOLDING
    ./analysis 4 -1 0 0 1 0 $ucf

    # echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf DATA PLOTTING TOGETHER
    # ./analysis 5 -1 0 0 1 0 $ucf    
done

rm $FNAME
