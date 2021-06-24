#!/bin/bash

for f in `find blif/*.blif`
do 
	echo $f
	./map -k 4 $f output.blif
	~/abc/abc -c "cec $f output.blif" # Equivalence check
	# echo -n '#Level: '
	# ~/abc/abc -c "read_blif output.blif;print_stats;" | grep 'lev =' | awk '{print $21}' # Level
	# echo -n '#LUT: '
	# grep .name output.blif | wc -l # LUT
done

# for f in aoi_9symml.blif 10aoi_alu4.blif 10aoi_big2.blif 10aoi_C1355.blif 10aoi_C6288.blif 10aoi_cht.blif 10aoi_cm138a.blif 10aoi_des.blif 10aoi_i2.blif 10aoi_i3.blif 10aoi_i4.blif 10aoi_k2.blif 10aoi_sample01.blif 10aoi_sample02.blif 10aoi_z4ml.blif
# do
#     echo $f
# done

