#|-lr_trimComplex.tcl :
#|  -generate simplified versions of ligand-protein complexes .
#|  -generates inputs for psfgen and namd .
#|  -based on the old script extAChEFrag.tcl .
#|  -author :-Carlos Z. Gómez Castro ;
#|  -date :
#|    -created :-2025-01-13.Mon ;
#|    -modified :-2025-01-14.Tue ;;
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
#|        -'src', 'source', 'input', 'inputSource' :
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
#|          -see also the 'source' arg .
#|          -'pdbid' :
#|            -specs common for both lig, and rec molecules ;
#|          -'pdbidL' :
#|            -specs for the ligand molecule ;
#|          -'pdbidR' :
#|            -specs for the receptor molecule ;
#|          -default value :
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
#|            -valid selection text for the VMD's command 'atomselect' ;
#|          -default value :-"" ;;
#|        -'selTxtR' :
#|          -specifies the VMD's selTxt for the ligand .
#|          -used optionally to restrict (exclude) receptor atoms interacting
#|           _ with the ligand .
#|          -acceptable values :
#|            -valid selection text for the VMD's command 'atomselect' ;
#|          -default value :-"" ;;
#|        -'atmSelL', 'atomSelLig' :
#|          -atom selection previously generated using the VMD's command
#|           _'atomselect' .
#|          -default value :-"" ;;
#|        -'atmSelR', 'atmSelCAR', 'atomSelRec' :
#|          -atom selection previously generated using the VMD's command
#|           _'atomselect' .
#|          -default value :-"" ;;
#|        -'cutoff', 'distCutoff', 'distance', 'dist' :
#|          -maximum distance between lig and rec atoms to be considered as
#|           _ interaction .
#|          -default value :-4.5 ;;
#|        -'gapMax', 'resGap', 'gap', 'gaps', 'gapSize', 'allowedGaps' :
#|          -max number of consecutive non-lig-interacting rec residues to
#|           _ be included in the final rec model .
#|          -gap residues connect more isolated rec fragments to yield bigger
#|           _ fragments .
#|          -acceptable values :-positive integers including zero .
#|          -recommended range of values :-0 to 5 ;
#|          -default value :-3 ;;
#|        -'tails', 'tailsLength', 'lengthTails', 'tailsLen', 'lenTails' :
#|          -max number of non-lig-interacting rec residues at the N- or
#|           C-terminations of each chain in the rec to be considered .
#|          -will aplly for each peptide chain in the rec .
#|          -peptide tails help to make more robust rec models .
#|          -acceptable values :-positive integers including zero .
#|          -recommended range of values :-0 to 5 ;
#|          -default value :-3 ;;
#|        -'tailN', 'NTail' :
#|          -same as 'tails' arg but specific for the peptide N-terminus ;
#|        -'tailC', 'CTail' :
#|          -same as 'tails' arg but specific for the peptide C-terminus ;
#|        -'mutateGaps', 'gapsMutate', 'mutate' :
#|          -specifies whether gap residues shold be mutated to Alanine .
#|          -glycine gap-rec residues will be allways ignored (not mutated) .
#|          -note: user-specified mutations may be considered in future
#|           _ versions .
#|          -acceptable values :-0 .-1 ;
#|          -default value :-1 ;
#|        -'keepProline', 'keepGapPro', 'proGapKeep' :- ;
#|        -'chargeTarget', 'targetCharge', 'mutateChargeTarget' :- ;
#|        - ;;
#|    -notes :
#|      - ;;
proc lr_trimComplex {complexType args} {
# global variables

# configuring logLib to manage log messages
  namespace import ::logLib::*
  set procName [lindex [info level 0] 0]
  set bll1 2; set bll2 [expr {$bll1 + 1}]   ;# base log level

# default value for variables
  set src "id"
  set pdbid "top"
  set pdbidL "top"
  set pdbidR "top"
  set prefix "auto"
  set selTxtL ""
  set selTxtR ""
  set atmSelL ""
  set atmSelR ""
  set cutoff 4.5
  set gapMax 3


  }   ;# proc lr_trimComplex

#|- ;




