#|-color_dlg.tcl :
#|  -collection of procedures to show ad4 results in VMD reps .
#|  -authors :
#|    -Carlos Z. Gómez Castro ;
#|  -date :
#|    -created :-2026-05-14.Thu ;
#|    -modified :-2026-05-14.Thu ;;
#|  -version :
#|    -001 :
#|      -procs color_lrcd, color_deltaG, color_drawFrames, and
#|       _ color_drawRestore recovered from backups (molInfo_tcrmaits.tcl) ;;;


proc color_lrcd {idl_ligs} {
  set min [[atomselect [lindex ${idl_ligs} 0] "index 0" frame 0] get user]
  set max $min
  foreach id ${idl_ligs} {
    set nf [molinfo $id get numframes]
    for {set f 0} {$f < $nf} {incr f} {
      set lrcdi [[atomselect $id "index 0" frame $f] get user]
      if {$lrcdi > $max} {set max $lrcdi}
      if {$lrcdi < $min} {set min $lrcdi}
      }
    }
  foreach id ${idl_ligs} {
    mol modcolor 1 $id User
    mol scaleminmax $id 1 $min $max
    }
  puts "min:max: $min:$max"
  }
# mol drawframes $id 1 "0:$nfmax"

proc color_deltaG {idl_ligs} {
  global molInfo selInfo
#  set min [lindex $selInfo($molInfo([lindex ${idl_ligs} 0],name),deltaG) 0]
  set min [lindex $molInfo([lindex ${idl_ligs} 0],deltaG) 0]
  set max $min
  foreach id ${idl_ligs} {
    set nf [molinfo $id get numframes]
    for {set f 0} {$f < $nf} {incr f} {
#      set delG [lindex $selInfo($molInfo($id,name),deltaG) $f]
      set delG [lindex $molInfo($id,deltaG) $f]
      [atomselect $id all frame $f] set user $delG
      if {$delG < $min} {set min $delG}
      if {$delG > $max} {set max $delG}
      }
    }
  foreach id ${idl_ligs} {
    mol modcolor 1 $id User
    mol scaleminmax $id 1 $min $max
    }
  puts "min:max: $min:$max"
  }

proc color_drawFrames {idl_ligs} {
  foreach id ${idl_ligs} {
    set nf [molinfo $id get numframes]
    mol drawframes $id 1 "0:$nf"
    mol on $id
    }
  }

proc color_drawRestore {idl_ligs} {
  foreach id ${idl_ligs} {
    mol modcolor 1 $id Name
    mol drawframes $id 1 "now"
    mol off $id
    }
  mol on [lindex ${idl_ligs} 0]
  mol top [lindex ${idl_ligs} 0]
#  display resetview
  }


#|- ;

