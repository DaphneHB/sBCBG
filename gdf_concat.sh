# To apply from the simulation's directory (log folder parent folder)
path=`pwd`
echo "Starting concatenation of GDF files for $path"

#if [ ! -f ./log/*0.gdf ]; then
#    echo "\tNo file to concatenate ---> Nothing changed"
#    exit
#fi


targetName="bigGDF.gdf"
echo -n "" > log/$targetName
echo "\tConcatenation in $targetName"
for N in "MSN" "FSI" "STN" "GPe" "GPi" ; do
    echo "====> $N data files concatenation" ;
    filename="$N-catFile.gdf"
    echo -n "" > $filename
    for f in `find log -type f -iname "$N-[0-9]*-[0-9].gdf"`; do
    	echo "\t....copying $f"
    	if [ "$f" != log/$filename ] ; then
    	    cat $f >> $filename
    	    rm -f $f
    	fi
    done
    echo "\n$N" >> log/$targetName
    sort -n $filename >> log/$targetName
    rm -f $filename
done

# getting rid of antag file if existing
for antag in "GPi_NMDA" "GPi_AMPA" "GPi_GABAA" "GPi_All" "GPi_NMDA+AMPA" "GPe_GABAA" "GPe_NMDA" "GPe_AMPA" "GPe_AMPA+GABAA" ; do
    for N in "MSN" "STN" "FSI" "GPe" "GPi" ; do
	prename="$antag""_$N"
	prefilename="$prename-catFile.gdf"
	echo -n "" > $prefilename
	echo "====> $N with $antag data files concatenation" ;
	for f in `find log -type f -iname "$prename-[0-9]*-[0-9].gdf"`; do
	    echo "\t....copying $f"
    	    if [ "$f" != log/$prefilename ] ; then
    		cat $f >> $prefilename
    		rm -f $f
    	    fi
	done
	echo "\n$prename" >> log/$targetName
	sort -n $prefilename >> log/$targetName
	rm -f $prefilename
    done
done

echo "DONE : big gdf file created & others gdf files removed"
