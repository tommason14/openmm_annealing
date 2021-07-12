set sel [atomselect top {not (name "DP.*" or name "VS.*")}]
$sel writepdb cleaned.pdb
exit
