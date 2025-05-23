#|-lrcd-ref-dlg.tcl :
#|  -declares procedures to calculate lrcd and other parameters associated 
#|   _ to docking calculations .
#|  -intended to create an output .csv file with energy, distCOM, and lrcd 
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

#|    -proc distCOMref {} :
#|      -calculate the distance between de center of mass of each conformation
#|       _ in an AD4 .dlg file relative and the COM of a ref compound .
#|      -arguments :
#|        -ref_selIds :- ;
#|        -l_selId :- ; ;;
proc distCOMref {ref_selIds l_selId {weights mass} {loSt stdout}} {
  
  }

#|  -namespace calc_lrcd :
#|    -collects different procedures to calculate lrcd and other ligand-receptor
#|     _ interaction properties ;;
namespace eval calc_lrcd {

#|    -global variables :
#|      -selInfo, molInfo ;
  global selInfo, molInfo

#|    -external libraries (within namespace) :
#|      -lrcd ;
  source lrcd.tcl

#|    -import :
#|      -::logLib::* ;
  namespace import ::logLib::*
  
#|    -export :
#|      - ;

#|    -variables :
#|      -l_refSelId :
#|        -list of selIds for the reference complexes .
#|        -for the moment will reference keys of the selInfo array ;
  variable l_refSelId {}
#|      -calc_lrcd_logFileName :- ;
  variable calc_lrcd_logFileName "log_calc_lrcd.txt"
#|      -calc_lrcd_outFileName :- ;
  variable calc_lrcd_outFileName "lrcd_summary.csv"

#|      -calc_lrcd_cutoff :- ;
  variable calc_lrcd_cutoff 4.5
#|      -calc_lrcd_components :- ;;
  variable calc_lrcd_components {"dct" "cmd" "lapa"}

#|    -namespace procs :
#|      - ;

  }   ;# namespace eval calc_lrcd

      - ;;
#|    -proc lrcd_summary_ref {} :
#|      -prints out a summary of docking results for a set of AD4 .dlg files .
#|      -considers an user-specified reference complex .
#|      -table of results (.csv format): :
#|        -header columns :
#|          -consecutive number (i) :- ;
#|          -current molecule id (id) :- ;
#|          -conformation frame (frm) :- ;
#|          -estimated free energy of binding (deltaG) :- ;
#|          -distance conformation COM to ref COM (distCOM) :- ;
#|          -ligand-receptor contact distance (lrcd) :- ;;
#|      -arguments :
#|        -ref_selIds :- ;
#|        -l_selId :- ;
#|        -args (optional arguments) :
#|          -lrcdType :- ;
#|          - ;;
#|      -notes :
#|        -requires the global array selInfo ;;;
proc lrcd_summary_ref {} {
  global selInfo
# loop over each reference
  foreach ref_selId $ref_selIds {
# get complex reference information
    set refId $selInfo($ref_selId,recSelId)
    
    }
  }


#|  - ;
#|- ;

