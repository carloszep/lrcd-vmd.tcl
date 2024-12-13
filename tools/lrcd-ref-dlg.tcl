#|-lrcd-ref-dlg.tcl :
#|  -declares procedures to calculate lrcd and other parameters associated 
#|   _ to docking calculations .
#|  -intended to create an output .csv file with energy, distCOM, and lrcd 
#|   _ for all the conformations in a group of .dlg files with a common ref .
#|  -requires library lrcd.tcl v. 1.1.3, and userInfo v. 0.0.6 .
#|  -requires preferentially a molInfo script already processed .
#|  -author :-Carlos Z. Gómez Castro ;
#|  -date :-2024-12-11.wed ;
#|  -version information :
#|    -version :-001 ;
#|    -changes in current version :
#|      -initial definitions as functions an procedures ;
#|    -previous changes :- ;;
#|  -procedures :


#|    -proc distCOMref {} :
#|      -calculate the distance between de center of mass of each conformation
#|       _ in an AD4 .dlg file relative and the COM of a ref compound .
#|      -arguments :
#|        -ref_selIds :- ;
#|        -l_selId :- ; ;;
proc distCOMref {ref_selIds l_selId {weights mass} {loSt stdout}} {
  
  }


#|    -proc lrcd_summary_ref {} :
#|      -prints a summary of docking results for a set of AD4 .dlg files .
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
proc lrcd_summary {} {
  }


#|  - ;
#|- ;
