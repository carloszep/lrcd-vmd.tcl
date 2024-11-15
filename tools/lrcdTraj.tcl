#|-lrcdTraj.tcl :
#|  -script to claculate LRCD along an MD trajectory ;
#|- ;

loadLib userInfo
loadLib anMD
loadLib lrcd

global selInfo trajInfo userInfo_version anMD_version
namespace import trajSelId::getFrame trajSelId::getDatX

puts "calculating LRCDi for a trajectory"
set ligSelTxt "resname LIG and noh"
set recSelTxt "protein and noh"
set logFilePref "log_lrcdTrajComp"
set trajStep 10

set id 0
set idRef 1

# output log file
set out [open "${logFilePref}_M3-M5.txt" w]

setSelId ligM3 molId $id selTxt $ligSelTxt step $trajStep updateSel 1 loSt $out
setSelId recM3 molId $id selTxt $recSelTxt step $trajStep updateSel 1 loSt $out
setSelId ligM5 molId $idRef selTxt $ligSelTxt step $trajStep updateSel 1 loSt $out
setSelId recM5 molId $idRef selTxt $recSelTxt step $trajStep updateSel 1 loSt $out

# collect selInfo for LRCD calculation
set selInfo(bclxl-M3,ligId) $id
set selInfo(bclxl-M3,ligSelTxt) $ligSelTxt
set selInfo(bclxl-M3,recId) $id
set selInfo(bclxl-M3,recSelTxt) $recSelTxt
set selInfo(bclxl-M3,title) "M3 in complex with BCL-XL"


# collect selInfo for LRCD calculation
set selInfo(bclxl-M5,ligId) $idRef
set selInfo(bclxl-M5,ligSelTxt) $ligSelTxt
set selInfo(bclxl-M5,recId) $idRef
set selInfo(bclxl-M5,recSelTxt) $recSelTxt
set selInfo(bclxl-M5,title) "M5 in complex with BCL-XL"

animate goto [simFrame 0 $idRef name]

# calculating lrcdMat array with all 6 components for reference complex
set lrcdMat {}
set lrcdMat [lrcdMat_dct $lrcdMat bclxl-M5 4.5 $out]
set lrcdMat [lrcdMat_cmd $lrcdMat bclxl-M5 4.5 $out]
set lrcdMat [lrcdMat_lapa $lrcdMat bclxl-M5 4.5 $out]
# calcualting lrcdVec array for t = 0 ns
set lrcdVecRef [lrcdVec_sum $lrcdMat $selInfo(bclxl-M5,title) $out]

# run over trajectory
set plots {}
set l_xDat {}
set l_yDat {}
mol top $id
trajSelId::init recM3
while {[trajSelId::iterate recM3 log $out]} {
  animate goto [getFrame]
# calculating lrcdMat array with all 6 components for each frame
  set lrcdMat {}
  set lrcdMat [lrcdMat_dct $lrcdMat bclxl-M3 4.5 "none"]
  set lrcdMat [lrcdMat_cmd $lrcdMat bclxl-M3 4.5 "none"]
  set lrcdMat [lrcdMat_lapa $lrcdMat bclxl-M3 4.5 "none"]
# calcualting lrcdVec array for current frame
  set lrcdVec [lrcdVec_sum $lrcdMat $selInfo(bclxl-M3,title) $out]

# calculate lrcd for the current frame
  set lrcdM3 [lrcnd $lrcdVecRef $lrcdVec "lrcdi for frame [getFrame]" "none"]
puts "time: [getDatX]; frame [getFrame]; lrcd $lrcdM3"

  lappend l_xDat [getDatX]
  lappend l_yDat $lrcdM3
  }

lappend plots [list $l_xDat $l_yDat "M3 - BCL-XL"]

# run over M5 sim
set l_xDat {}
set l_yDat {}
mol top $idRef
trajSelId::init recM5
while {[trajSelId::iterate recM5 log $out]} {
  animate goto [getFrame]
# calculating lrcdMat array with all 6 components for each frame
  set lrcdMat {}
  set lrcdMat [lrcdMat_dct $lrcdMat bclxl-M5 4.5 "none"]
  set lrcdMat [lrcdMat_cmd $lrcdMat bclxl-M5 4.5 "none"]
  set lrcdMat [lrcdMat_lapa $lrcdMat bclxl-M5 4.5 "none"]
# calcualting lrcdVec array for current frame
  set lrcdVec [lrcdVec_sum $lrcdMat $selInfo(bclxl-M5,title) $out]

# calculate lrcd for the current frame
  set lrcdM5 [lrcnd $lrcdVecRef $lrcdVec "lrcdi for frame [getFrame]" "none"]
puts "time: [getDatX]; frame [getFrame]; lrcd $lrcdM5"

  lappend l_xDat [getDatX]
  lappend l_yDat $lrcdM5
  }

lappend plots [list $l_xDat $l_yDat "M5 - BCL-XL"]

# creating xmgrace plot
agrXY $plots "lrcdTrajComp" \
  title "Trajectory LRCND" \
  xLabel "Time (ns)" yLabel "LRCND" loSt $out

close $out

