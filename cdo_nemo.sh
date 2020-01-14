##################################################################################################################################
## This bash script searches in all directories and subdirectories of a specified path monthly NEMO model outputs (each month is
## stored in a different directory). The script then uses CDO to select a group of variables and regrid them from the 
## ocean model native grid to a regular lat-lon grid. 
## Note that the script has been designed to be run (also) on incomplete/running simulations. If an alredy remapped file is found 
## in any directory, after making sure the file is the right size (i.e. it is not corrupted) the scripts will skip that file and
## move on to the next directory until it finds new unprocessed model outputs. 
##################################################################################################################################

#!/bin/bash

for file in `find /gws/nopw/j04/pmip4_vol1/users/vittoria/u-bk453/ -maxdepth 2 -type f -name 'nemo_bk453o_1m_*-T.nc'`;
 do 
    ##get the path to the directory for each file read in
    datadir=$(dirname $file);
    ##check if the file is alredy there (-e 'exists')) from previous runs of the code (and if it has the right size!)
   if [ -e $datadir/ThetaoSo_$(basename $file) ]; then
    filesize=$(wc -c < $datadir/ThetaoSo_$(basename $file) ) 
   else
   filesize=0 ;fi
 
    if [ $filesize -ge 58389737 ]; then
    echo 'filesize=' $filesize
    echo 'Found File - skipping '$file'';
    else
    echo 'filesize=' $filesize
    echo 'processing '$file''
    ##submit one separate job for each file found in $datadir. The output file is stored in $datadir and it has same name as input 
    ##file but with a prefix "prefix_"
     bsub -W 00:20 -J 'cdo' cdo -L -remap,global_lonlat_1deg,remap_Tgrid_weigths.nc  -selname,so,thetao,thkcello   $file $datadir/ThetaoSo_$(basename $file)
     
   fi  
 done 
