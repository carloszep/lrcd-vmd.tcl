#|-lr_trimComplex.tcl :
#|  -generate simplified versions of ligand-protein complexes .
#|  -generates inputs for psfgen and namd .
#|  -based on the old script extAChEFrag.tcl .
#|  -author :-Carlos Z. Gómez Castro ;
#|  -date :
#|    -created :-2025-01-13.Mon ;
#|    -modified :-2025-01-13.Mon ;;
#|  -version :
#|    -001 :
#|      -procedure first defined ;;
#|  -source (external files) :
#|    -source logLib.tcl ;;
  source logLib.tcl

#|  -proc lr_trimComplex { } :
#|    -writes pdb files for each fragment in the protein that interacts with
#|     _ the ligand .
#|    -creates lists of residues to be mutated .
#|    -import namespaces :
#|      -::logLib::* ; ;
#|    -arguments :
#|      -'complexType' :
#|        -specifies the type of structure to be considered .
#|        -acceptable values :
#|          -"pdb" :
#|            -both the ligand and the receptor are parts of a PDB model .
#|            -the 'source', 'pdbid' arguments must be specified .
#|            -args 'pdbidL' and 'pdbidR' will be ignored ;
#|          -"dlg" :
#|            -the ligand is an AD4 .dlg file .
#|            -the receptor is in a different molecule .
#|            -*** this option is not implemented yet *** ;
#|          -"other" :
#|            -*** this option is not implemented yet *** ;;;
#|      -args (variable arguments) :
#|        -'source', 'input', 'inputSource' :
#|          -type of files to be loaded .
#|          -acceptable values :
#|            -"download", "pdbid", "pdbload" :
#|              -pdb file(s) downloaded using the VMD command 'pdbload' .
#|              -arg 'pdbid' will be interpreted as a 4-letter PDBID ;
#|            -"file", "loadFile" :
#|              -file names to locally load the files are assumed .
#|              -can contain absolute or relative file paths ;
#|            -"id", "ids", "molId", "molIds", "loaded", "vmdId" :
#|              -the molecules are already loaded into VMD .
#|              -the VMD id(s) should be specified in args 'pdbid', 'pdbidl',
#|               _ or 'pdbidr' ;;
#|          -default value :
#|            -"id" ;;
#|        -'pdbid', 'pdbidL', 'pdbidR' :
#|          -can be interpreted as a 4-letter PDBID, a file name, or a VMD Id,
#|           _ depending on the 'source' argument .
#|          -specs for a common, ligand, and receptor molecules, resp. .
#|          -see also the 'source' arg .
#|          -acceptable values :
#|            -'top' :-use to specify the VMD's top molecule ;;;
#|        -'prefix', 'outPrefix', 'outputPrefix' :
#|          -prefix used for all pdb file generated .
#|          -the specified name may include abs or rel file paths .
#|          -acceptable values :
#|            -any string acceptable as part of a file name .
#|            -"auto" :- ;;
#|          -default value :-"auto" ;;
#|        -'selTxtL' :
#|          -specifies the VMD's selTxt for the ligand .
#|          -acceptable values :
#|            - ;
#|          -default value :- ;;
#|        -'selTxtR' :
#|          -specifies the VMD's selTxt for the ligand .
#|          -used optionally to restrict (exclude) receptor atoms interacting
#|           _ with the ligand .
#|          -acceptable values :
#|            - ;;
#|        -'atmSelL', 'atomSelLig' :
#|          - ;
#|        -'atmSelR', 'atmSelCAR', 'atomSelRec' :- ;
#|        -'cutoff', 'distCutoff', 'distance', 'dist' :- ;
#|        -'gapMax', 'resGap', 'gap', 'gaps', 'gapSize', 'allowedGaps' :- ;
#|        -'tailsLength', 'lengthTails', 'tailsLen', 'lenTails' :- ;
#|        -'tailN', 'NTail' :- ;
#|        -'tailC', 'CTail' :- ;
#|        -'mutateGaps', 'gapsMutate', 'gapsMutateAla', 'mutate' :- ;
#|        -'keepProline', 'keepGapPro', 'proGapKeep' :- ;
#|        -'chargeTarget', 'targetCharge', 'mutateChargeTarget' :- ;
#|        - ;;
#|    -notes :
#|      - ;;
proc lr_trimComplex {complexType  } {
# global variables

# configuring logLib to manage log messages
  namespace import ::logLib::*
  set procName [lindex [info level 0] 0]
  set bll1 2; set bll2 [expr {$bll1 + 1}]   ;# base log level

# default value for variables
  


  }   ;# proc lr_trimComplex

#|- ;

