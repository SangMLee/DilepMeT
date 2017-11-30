import DMsimp_s_spin1-vector
generate p p > xd xd~ l+ l- [QCD]
add process p p > xd xd~ l+ l- j [QCD]
output DMsimp_s_spin1-vector_1j_qcd
launch
shower=ON
done
./Cards/param_card_DMs1.dat
./Cards/run_card_DMs1.dat
./Cards/shower_card_DMs1_1j.dat
done
# Will need to run delphes manually for the loop guys
# zcat <hepmc> | singularity exec ~/Images/CCMadgraph.img /code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC <config> <output> -

import DMsimp_s_spin1-vector
generate p p > xd xd~ l+ l- [QCD]
output DMsimp_s_spin1-vector_0j_qcd
launch
shower=ON
done
./Cards/param_card_DMs1.dat
./Cards/run_card_DMs1.dat
./Cards/shower_card_DMs1_0j.dat
done
# Match the output name of the regular guys
# zcat DMsimp_s_spin1-vector_0j_qcd/Events/run_01/events_PYTHIA8_0.hepmc.gz | singularity exec ~/Images/CCMadgraph.img /code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC ./DMsimp_s_spin1-axial_0j/Cards/delphes_card_CMS.dat DMsimp_s_spin1-vector_0j_qcd/Events/run_01/tag_1_delphes_events.root -

import DMsimp_s_spin1-vector
generate p p > xd xd~ l+ l- [QCD]
add process p p > xd xd~ l+ l- j [QCD]
output DMsimp_s_spin1-vector_1j_qcd_Qcut100
launch
shower=ON
done
./Cards/param_card_DMs1.dat
./Cards/run_card_DMs1.dat
./Cards/shower_card_DMs1_1j_Qcut100.dat
done

# Match the output name of the regular guys
# zcat DMsimp_s_spin1-vector_0j_qcd/Events/run_01/events_PYTHIA8_0.hepmc.gz | singularity exec ~/Images/CCMadgraph.img /code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC ./DMsimp_s_spin1-axial_0j/Cards/delphes_card_CMS.dat DMsimp_s_spin1-vector_0j_qcd/Events/run_01/tag_1_delphes_events.root -

# import DMsimp_s_spin1-vector
# generate p p > xd xd~ l+ l- [QCD]
# add process p p > xd xd~ l+ l- j [QCD]
# output DMsimp_s_spin1-vector
# launch
# shower=ON
# done
# ./CARDS/run_card.dat
# ./CARDS/shower_card.dat
# ./CARDS/param_card.dat
# done

# 0j, no shower/matching changes needed
import DMsimp_s_spin1-vector
generate p p > xd xd~ l+ l-
output DMsimp_s_spin1-vector_0j
launch
shower=PYTHIA8
delphes=ON
done
./Cards/param_card_DMs1.dat
done


# 0j, no shower/matching changes needed
import DMsimp_s_spin1-vector
generate p p > xd xd~ e+ e-
output DMsimp_s_spin1-vector_ee_0j
launch
shower=PYTHIA8
delphes=ON
done
./Cards/param_card_DMs1.dat
done

# 0j, no shower/matching changes needed
import DMsimp_s_spin1-vector
generate p p > xd xd~ mu+ mu-
output DMsimp_s_spin1-vector_mm_0j
launch
shower=PYTHIA8
delphes=ON
done
./Cards/param_card_DMs1.dat
done

# 0j, no shower/matching changes needed
import DMsimp_s_spin1-axial
generate p p > xd xd~ l+ l-
output DMsimp_s_spin1-axial_0j
launch
shower=PYTHIA8
delphes=ON
done
./Cards/param_card_DMs1axial.dat
done

# 0j, no shower/matching changes needed
import DMsimp_s_spin1-axial
generate p p > xd xd~ mu+ mu-
output DMsimp_s_spin1-axial_mm_0j
launch
shower=PYTHIA8
delphes=ON
done
./Cards/param_card_DMs1axial.dat
done


# 0j, no shower/matching changes needed
import DMsimp_s_spin1-axial
generate p p > xd xd~ e+ e-
output DMsimp_s_spin1-axial_ee_0j
launch
shower=PYTHIA8
delphes=ON
done
./Cards/param_card_DMs1axial.dat
done

# 1j try default matching
import DMsimp_s_spin1-vector
generate p p > xd xd~ l+ l-
add process p p > xd xd~ l+ l- j
output DMsimp_s_spin1-vector_1j
launch
shower=PYTHIA8
delphes=ON
done
./Cards/param_card_DMs1.dat
done

