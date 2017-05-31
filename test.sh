for antag in "GPi_NMDA" "GPi_AMPA" "GPi_GABAA" "GPi_All" "GPi_NMDA+AMPA" "GPe_GABAA" "GPe_NMDA" "GPe_AMPA" "GPe_AMPA+GABAA" ; do
    echo "====> $antag data files concatenation" ;
    filename="$antag-catFile.gdf"
    echo "filename = $filename"
    for N in "MSN" "STN" "FSI" "GPe" "GPi" ; do
	prename="$antag""_$N"
	echo "${prename}" 
	for f in `find log -type f -iname "$prename-[0-9]*-[0-9].gdf"`; do
	    echo "the file : $f"
	done
    done
done

