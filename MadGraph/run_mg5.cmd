import DMsimp_s_spin1-vector
generate p p > xd xd~ l+ l- [QCD]
add process p p > xd xd~ l+ l- j [QCD]
output DMsimp_s_spin1-vector
launch
shower=ON
done
./CARDS/run_card.dat
./CARDS/shower_card.dat
./CARDS/param_card.dat
done

