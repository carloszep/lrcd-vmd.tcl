#|-lrcd-ref-dlg.tcl :
#|  -declares procedures to calculate lrcd and other parameters associated 
#|   _ to docking calculations .
#|  -is intended to create an output .csv file with energy, distCOM, and lrcd 
#|   _ for all the conformations in a group of .dlg files with a common ref .
#|  -requires library lrcd v. 1.1.3, and userInfo v. 0.0.6 .
#|  -requires preferentially a molInfo script already processed .
#|  -author :-Carlos Z. Gómez Castro ;
#|  -date :-2024-12-11.wed ;
#|  -version :
#|    -001 :
#|      -initial definitions as functions an procedures ;
#|  -notes :
#|    - ;

#|  -external libraries :
#|    -logLib .
#|    -userInfo .
#|    -lrcd ;
source logLib.tcl
source userInfo.tcl

#|  -namespace calc_lrcd :
#|    -collects different procedures to calculate lrcd and other ligand-receptor
#|     _ interaction properties ;;
namespace eval calc_lrcd {

#|    -external libraries (within namespace) :
#|      -lrcd ;
  source lrcd.tcl

#|    -import :
#|      -::logLib::* ;
  
#|    -export :
#|      - ;

#|    -variables :
#|      -l_refSelId :
#|        -list of selIds for the reference complexes .
#|        -for the moment will reference keys of the selInfo array ;
  variable l_refSelId {}
#|      -calc_lrcd_logFileName :- ;
  variable calc_lrcd_logFileName "log_calc_lrcd.txt"

#|      -calc_lrcd_cutoff :- ;
  variable calc_lrcd_cutoff 4.5
#|      -calc_lrcd_components :- ;;
  variable calc_lrcd_components {"dct" "cmd" "lapa"}

#|    -namespace procs :
#|      - ;
  proc resultsTable {} {
    }

  }   ;# namespace eval calc_lrcd

      - ;;

#|- ;
