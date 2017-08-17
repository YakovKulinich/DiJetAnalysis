#/bin/bash!
FNAME=uc.txt

# uncertainties to run over
# UC=( -3 -2 0 2 3 )#
#UC=( -1 -2 -3 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 0 1 2 3 9 10 11 12 13 14 15 16 17 18 )
#UC=( -1 -2 -3 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 0 16 17 )
UC=( -9 -10 -11 -12 -13 -14 -16 -19 0 9 13 14 16 19 20 )
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
    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf MC JETS
 #   ./analysis 0 -1 0 0 0 0 $ucf
    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf MC HISTOS
 #   ./analysis 1 -1 0 0 0 0 $ucf
    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf DATA UNFOLDING
    ./analysis 2 -1 0 0 1 0 $ucf
    echo !!!!!!!!!!!!!!!!!!!!!!!!!!!! $ucf DATA PLOTTING TOGETHER
    ./analysis 4 -1 0 0 1 0 $ucf    
done

rm $FNAME
