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
    	echo "\t....analyzing $f"
    	if [ "$f" != log/$filename ] ; then
    	    cat $f >> $filename
    	    rm -f $f
    	fi
    done
    echo "$N" >> log/$targetName
    sort -n $filename >> log/$targetName
    rm -f $filename
done
echo "DONE : big gdf file created & others gdf files removed"
