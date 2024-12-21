filename="offline_PinnedH2O.out" # FOM
#filename="rom39_PinnedH2O.out" # compare MD
#filename="39_force_PinnedH2O.out" # compare force

# Extracting H1, H2, F1, F2 from MGmgol output log
# if FOM, these files contain the FOM results
# if compare MD, these files contain the results with projected orbitals
awk '/H1 / {print $3, $4, $5}' $filename > H1_$filename
awk '/H2 / {print $3, $4, $5}' $filename > H2_$filename
awk '/F1 / {print $6, $7, $8}' $filename > F1_$filename
awk '/F2 / {print $6, $7, $8}' $filename > F2_$filename

# if compare force, files with "_fom" contain the FOM results
# files with "_rom" contain the results with projected orbitals
if [[ "$filename" == *"force_"* ]]; then
    sed -n '1~2p' H1_$filename > H1_rom$filename
    sed -n '1~2p' H2_$filename > H2_rom$filename
    sed -n '1~2p' F1_$filename > F1_rom$filename
    sed -n '1~2p' F2_$filename > F2_rom$filename

    sed -n '2~2p' H1_$filename > H1_fom$filename
    sed -n '2~2p' H2_$filename > H2_fom$filename
    sed -n '2~2p' F1_$filename > F1_fom$filename
    sed -n '2~2p' F2_$filename > F2_fom$filename
fi

rm -rf snapshot0_*
