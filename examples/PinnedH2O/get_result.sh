#filename="offline_PinnedH2O.out" # FOM
#filename="rom39_PinnedH2O.out" # ROM compare MD
#filename="39_force_PinnedH2O.out" # ROM compare force

#filename="PinnedH2O_ref.out" # FOM
filename="PinnedH2O_rom_3DOF_2_2_34.out" # ROM PinnedH2O 3DOF MD

# Extracting H1, H2, F1, F2 from MGmgol output log
# if FOM, these files contain the FOM results
# if compare MD, these files contain the results with projected orbitals
awk '/O1 / {print $4, $5, $6}' $filename > O1_$filename
awk '/H1 / {print $3, $4, $5}' $filename > H1_$filename
awk '/H2 / {print $3, $4, $5}' $filename > H2_$filename
awk '/O1 / {print $7, $8, $9}' $filename > f_O1_$filename
awk '/H1 / {print $6, $7, $8}' $filename > f_H1_$filename
awk '/H2 / {print $6, $7, $8}' $filename > f_H2_$filename

# if compare force, files with "_fom" contain the FOM results
# files with "_rom" contain the results with projected orbitals
if [[ "$filename" == *"force_"* ]]; then
    sed -n '1~2p' H1_$filename > H1_rom$filename
    sed -n '1~2p' H2_$filename > H2_rom$filename
    sed -n '1~2p' f_H1_$filename > f_H1_rom$filename
    sed -n '1~2p' f_H2_$filename > f_H2_rom$filename

    sed -n '2~2p' H1_$filename > H1_fom$filename
    sed -n '2~2p' H2_$filename > H2_fom$filename
    sed -n '2~2p' f_H1_$filename > f_H1_fom$filename
    sed -n '2~2p' f_H2_$filename > f_H2_fom$filename
fi

rm -rf snapshot0_*
