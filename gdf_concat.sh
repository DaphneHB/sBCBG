# To apply from the simulation's directory
path=`pwd`
echo "Starting concatenation of GDF files for $path"

if [ ! -f ./log/*0.gdf ]; then
    echo "\tNo file to concatenate ---> Nothing changed"
    exit
fi

targetName="bigGDF.gdf"
echo -n "" > $targetName
echo "\tConcatenation in $targetName"
for f in `find log -type f -iname '*.gdf'` ; do
    if [ "$f" != log/$targetName ] ; then
	cat $f >> $targetName
	rm -f $f
    fi
done
sort -n $targetName > log/$targetName
rm -f $targetName
echo "\tDONE : big gdf file created & others gdf files removed"
