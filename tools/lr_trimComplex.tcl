#|-lr_trimComplex.tcl :
#|  -generate simplified versions of ligand-protein complexes .
#|  -generates inputs for psfgen .
#|  -based on the old script extAChEFrag.tcl .
#|  -author :-Carlos Z. Gómez Castro ;
#|  -public repository :-https://github.com/carloszep/lrcd-vmd.tcl ;
#|  -date :
#|    -created :-2025-01-13.Mon ;
#|    -modified :-2025-03-25.Thu ;;
set lr_trimComplex_version 008
#|  -version :
#|    -008 :
#|      -variable argument 'regenerateCad' replaced by 'userPatchesCad' .
#|      -bug fixed for argument userPatchesCad ;
#|      -no tested yet ;
#|    -007 :
#|      -PRO residues (besides GLY) are always excluded from the QM region .
#|      -all arguments of apply_restrictions have default value .
#|      -added optional argument pgnOutPath .
#|      -redirected some output ;
#|    -006 :
#|      -proc apply_restrictions added and tested .
#|      -routine for loading the created structure (after psfgen) to: :
#|        -write structures with restrictions for equil, MD sims, and QM/MM ;
#|      -for the moment only functionality for QM/MM ;
#|    -005 :
#|      -calling command reported in the header of the pgn script .
#|      -incorporated selTxtL_noInt argument .
#|      -non-interacting ligands can be added to the complex without growing-up
#|       _ receptor model residues due to interacting ligands .
#|      -ligSel and interactCASel eliminated as arguments  ;
#|    -004 :
#|      -bugfix for the first segment's tail .
#|      -multiple ligand selections are acceptable as a list .
#|      -definition of ligand segments implemented .
#|      -writting pdb model in file baseLR.tmp for psfgen to read coordinates ;
#|    -003 :
#|      -first working version to define segments .
#|      -rewritten as namespace ;
#|    -002 :
#|      -version almost completed but missing tails .
#|      -implementing new code for variable gap lengths ;
#|    -001 :
#|      -procedure's command-line interface .
#|      -procedure first defined .
#|      -original code tested ;
#|    -to do :
#|      -to join segments contiguous or overlapping due to tails .
#|      -to fix accumulation of string variables on repeated use
#|       _ (reinitialization) .
#|      -to implement the chargeTarget function .
#|      -to implement managing lig and rec molecules from different files .
#|      -manage crystallographic water .
#|      -to fix direct use of psfgen package (option pgnWrite 0) .
#|      -add namespace procs to optionally pre-process input arguments .
#|      -refine namespace procedures ;;
#|  -source (external files) :
#|    -logLib.tcl :
#|      -the loadLib proc is defined in the .vmdrc file ;;;
loadLib logLib

#|-namespace eval lr_trimComplex :
#|  -creates psfgen structures with fragments of ligand-receptor complexes .
#|  -detects protein fragments that interacts with ligands .
#|  -modifies the resulted "trimed" complex by adding gaps, tails,
#|   _ and/or mutations  .
namespace eval lr_trimComplex {

#|  -import :
#|    -logLib::* ;
  namespace import ::logLib::*

#|  -export :
#|    -trimComplex ;
  namespace export lr_trimComplex apply_restrictions

#|  -arguments (lr_trimComplex proc) :
#|    -'complexType' :
#|      -specifies the type of structure to be processed .
#|      -acceptable values :
#|        -"pdb", "singleMol", "singlePDB" :
#|          -both the ligand and the receptor are parts of a PDB model .
#|          -'source' and 'pdbid' must be specified .
#|          -args 'pdbidL' and 'pdbidR' will be ignored ;
#|        -"selInfo", "selId" :
#|          -lig and rec information is taken form the global selInfo array .
#|          -the 'selId' must be specified .
#|          -the complex may come from a single molecule or from two ;
#|        -"dlg" :
#|          -the ligand is an AD4 .dlg file .
#|          -the receptor is in a different molecule .
#|          -*** this option is not implemented yet *** ;

#|  -variables (optional argments) :
#|    -'src', 'source', 'input', 'inputSource' :
#|      -type of files to be loaded .
#|      -acceptable values :
#|        -"download", "pdbId", "pdbLoad", "PDB" :
#|          -pdb file(s) downloaded using the VMD command 'pdbload' .
#|          -arg 'pdbid' will be interpreted as a 4-letter PDBID ;
#|        -"file", "loadFile", "localFile", "localPDB" :
#|          -file names to locally load the files are assumed .
#|          -can contain absolute or relative file paths ;
#|        -"id", "ids", "molId", "molIds", "loaded", "vmdId", "vmd" :
#|          -the molecules are already loaded into VMD .
#|          -the VMD id(s) should be specified in args 'pdbid', 'pdbidl',
#|           _ or 'pdbidr' ;;
#|      -default value :-"id" ;;
  variable src "id"

#|    -'selId', 'selInfo', 'selInfoId', 'selInfoKey' :
#|      -selId key to use the information stored in the selInfo array .
#|      -implies that the lig and rec molecules were
#|       _ already loaded into VMD (maybe through a molInfo file) .
#|      -the atom selections registered in the selInfo array are expected
#|       _ to be more general (i.e., "protein and chain A", "not protein") .
#|      -the atomselections are processed to detect interacting residues .
#|      -selInfo array keys considered :
#|        -'<selId>,ligId', '<selId>,recId',
#|         _ '<selId>,ligSelTxt', and '<selId>,recSelTxt' .
#|        -'<selId>,frame' (optionally) ;
#|      -variables to be assigned (respectively): pdbIdL, pdbIdR,
#|       _ selTxtL, and selTxtR .
#|      -arguments that could override selInfo-stored values: 'pdbId*',
#|       _ 'frame', 'selTxt*' .
#|      -acceptable values :
#|        -an already registered selId for the global selInfo array .
#|        -"" :-the 'selId' arg will be ignored ;;
#|      -default value :-"" ;;  
#|    -'frame', 'frm', 'index' :
#|      -for .dlg files, it allows to specify the desired conformation .
#|      -default value ;-0 ;;
  variable selId ""
  variable frame 0

#|    -'pdbId', 'pdbFile', 'id', 'vmdId', 'molId' :
#|      -specs common for both lig and rec molecules .
#|      -can be interpreted as a 4-letter PDBID, a file name, or a VMD Id,
#|       _ depending on the 'source' argument .
#|      -this arg along with 'pdbIdL' and 'pdbIdR' may override the
#|       _ specs gotten from the selInfo array .
#|      -default value :-"top" ;;
#|    -'pdbIdL', 'pdbFileL', 'idL', 'vmdIdL', 'molIdL' :
#|      -same as arg 'pdbId' but only applies for the ligand molecule .
#|      -default value :-"" ;;
#|    -'pdbIdR', 'pdbFileR', 'idR', 'vmdIdR', 'molIdR' :
#|      -same as arg 'pdbId' but only applies for the receptor molecule .
#|      -default value :-"" ;;
  variable pdbId "top"
  variable pdbIdL ""
  variable pdbIdR ""

#|    -'selTxtL', 'l_selTxtL', 'selTxtLs' :
#|      -specifies a list of VMD's atom selection text for the ligand(s) .
#|      -each sel text must correspond to a single ligand .
#|      -for a single sel text use always quotes and curly braces, i.e.
#|       _ '{"resname LIG"}' .
#|      -acceptable values :
#|        -a list of valid selection text for the VMD's command 'atomselect' ;
#|      -default value :-{"not (protein or nucleic or water)"} ;;
#|    -'selTxtR' :
#|      -specifies the VMD's atom selection text for the receptor atoms
#|      -overrides the atom selection text from the 'selInfo' array .
#|      -acceptable values :
#|        -valid selection text for the VMD's command 'atomselect' ;
#|      -default value :-"protein" ;;
#|    -'selTxtL_noInt', 'l_selTxtL_noInt', 'selTxtL_nonInteract',
#|     _ 'selTxt_nonInteractL', 'extraLig', 'extraLigands' :
#|      -list of VMD's atom selection texts for non-interacting ligands .
#|      -these ligands are not considered to detect receptor's
#|       _ ligand-interacting residues .
#|      -A segment is added for each ligand .
#|      -for a single sel text use always quotes and curly braces, i.e.
#|       _ '{"resname LIG"}' .
#|      -accepatable values :
#|        -"" :-no non-interacting ligands are added ;
#|        -a list of valid selection text for the VMD's command 'atomselect' ;
#|      -default value :-"" ;;
  variable selTxtL {"not (protein or nucleic or water)"}
  variable selTxtR "protein or nucleic"
  variable selTxtL_noInt ""

#|    -'ligSel', 'atmSelL', 'atomSelLig' :
#|      -atom selection pregenerated using 'atomselect' VMD's command .
#|      -overrides 'selInfo(<selId>,ligSelTxt)' or 'selTxtL' .
#|      -default value :-"" ;
#|      -note: discontinuing variable  ;
#|    -'interactCASel', 'iteractingCA', 'atmSelR', 'recSel',
#|     _ 'atomSelRec', 'atmSelCA' :
#|      -atom selection pregenerated using 'atomselect' VMD's command .
#|      -it must include only the alpha carbon atoms from the residues
#|       _ in the receptor in direct contact with the ligand .
#|      -overrides 'selInfo(<selId>,recSelTxt)' and 'selTxtR' .
#|      -default value :-"" ;
#|      -note: discontinuing variable  ;
#  variable ligSel ""
#  variable interactCASel ""

#|    -'l_topFile', 'l_topFiles', 'topFile', 'topFiles',
#|     _ 'l_topologyFile', 'topologyFiles' :
#|      -list of topology files to be used by psfgen .
#|      -the file names may include either relative or absolute paths .
#|      -this list will be read before the top files found in the 'topDir' .
#|      -default value :-{} ;;
#|    -'topDir', 'topPath', 'topologyDir', 'topologyPath' :
#|      -path to the directory where topology files can be found .
#|      -acceptable values :
#|        -an absolute or relative dir path (ending with "/") .
#|        -a list of dir paths is also acceptable .
#|        -"" :-no directory is considered to search topology paths ;;
#|      -default value :-"" ;;
#|    -'l_topExt', 'l_topFileExt' :
#|      -list of known file extensions for CHARMM topology .
#|      -default value :-{top rtf str} ;;
#|    -'pdbAliasCad' :
#|      -string with several lines for "pdbalias" psfgen commands .
#|      -default value :
#|        -"pdbalias residue HIS HSD\npdbalias atom ILE CD1 CD\npdbalias
#|         _ residue HOH TP3M\npdbalias atom TP3M O OH2\n" ;;
#|    -'userPatchesCad', 'usePatchCad', 'patchCad', 'patchesCad', 'patches',
#|     _ 'userPatches' :
#|      -string with one or more lines with user-specified patches .
#|      -see also psfgen user's manual .
#|      -for some partches will be neccesary to add the line
#|       _ "regenerate angles dihedrals" to correctly apply the patch(es) .
#|      -example: "patch ASPP P:209; regenerate angles dihdrals" ;
#|      -default value :-"" ;;
  variable l_topFile {}
  variable topDir ""
  variable l_topExt {top rtf str}
  variable pdbAliasCad "pdbalias residue HIS HSD\npdbalias atom ILE CD1 CD\npdbalias residue HOH TP3M\npdbalias atom TP3M O OH2\n"
  variable userPatchesCad ""

#|    -'prefix', 'outPrefix', 'outputPrefix' :
#|      -prefix used for all pdb file generated .
#|      -the specified name may include abs or rel file paths .
#|      -acceptable values :
#|        -any string acceptable as part of a file name .
#|        -"", "auto" :- ;;
#|      -default value :-"" ;;
#|    -'workPath', 'workFolder', 'workDir',
#|     _ 'outPath', 'outFolder', 'outDir' :
#|      -directory path to be prepended to output files (before prefix) .
#|      -the path will be created in case it is not found ;
#|      -acceptable values :
#|        -a valid either relative or absolute directory path .
#|        -"", ".", "./" :
#|          -the current folder is used as workPath ;;
#|      -defautl value :-"" ;
#|    -'pgnWrite', 'pgnScript', 'pgnFile', 'psfgenFile', 'psfgenScript' :
#|      -controls whether a psfgen is generated 'on-the-fly' or a psfgen
#|       _ config pgn file is written for later processing .
#|      -acceptable values :
#|        -0 :-the psfgen structures are created within this proc ;
#|        -1 :-a psfgen config pgn file is created instead ;;
#|      -default value :-1 ;;
#|    -'segPrefix', 'segmentPrefix', 'segLetter' :
#|      -one-to-three-letter name prefix for the created segments .
#|      -default value :-"F" ;;
#|    -'pgnOutPath', 'pgnDirOut', 'pgnDirOutput' :
#|      -path specified for the psfgen output .
#|      -where the .psf and .pdb files will be found after running the script
#|       _ generated by lr_trimComplex using psfgen .
#|      -default value :
#|        -"../struct/" :
#|        -it is supposed that lr_trimComplex is run in a psfgen/ dis, where
#|         _ the .pgn.tcl script will be run but the psfgen output will be
#|         _ at the struct/ dir at the same level as psfgen/ ;;;
  variable prefix ""
  variable workPath ""
  variable pgnWrite 1
  variable segPrefix "F"
  variable pgnOutPath "../struct/"

#|    -'cutoff', 'distCutoff', 'distance', 'dist', 'r' :
#|      -maximum distance between lig and rec atoms to be considered as
#|       _ interaction .
#|      -default value :-4.5 ;;
  variable cutoff 4.5

#|    -'gapMax', 'maxGap', 'resGap', 'gap', 'gaps', 'gapSize', 'gapLen' :
#|      -max number of consecutive non-lig-interacting rec residues to
#|       _ be included in the final rec model .
#|      -gap residues connect isolated rec frags to yield bigger fragments .
#|      -acceptable values :-positive integers including zero .
#|      -recommended range of values :-0 to 5 ;
#|      -default value :-3 ;;
#|    -'mutateGaps', 'gapsMutate', 'mutate' :
#|      -specifies whether gap residues should be mutated to residue 'gapMut' .
#|      -apply to tail residues too .
#|      -note: user-specified mutations may be considered in future versions .
#|      -acceptable values :-0 .-1 ;
#|      -default value :-1 ;
#|    -'gapMut', 'gapMutation', 'gapRes', 'gapResidue', gapResName :
#|      -specify which residue are gap residues mutated to .
#|      -acceptable values :-CHARMM resname ;
#|      -default value :-"ALA" ;;
  variable gapMax 3
  variable mutateGaps 1
  variable gapMut "ALA"

#|    -'tails', 'tailsLength', 'lengthTails', 'tailsLen', 'lenTails' :
#|      -max number of non-lig-interacting rec residues at the N- or
#|       C-terminations of each chain in the rec to be considered .
#|      -will aplly for each peptide chain in the rec .
#|      -peptide tails help to make more robust rec models .
#|      -acceptable values :-positive integers including zero .
#|      -recommended range of values :-0 to 5 ;
#|      -default value :-2 ;;
#|    -'tailN', 'NTail', 'tail-N', 'N-tail' :
#|      -same as 'tails' arg but specific for the peptide N-terminus ;
#|    -'tailC', 'CTail', 'tail-C', 'C-tail' :
#|      -same as 'tails' arg but specific for the peptide C-terminus ;
#|    -'internalTails', 'tailsInternal', 'intTails' :
#|      -specify whether tails are also added to segments within a chain
#|       _ or only at the extremes of each chain .
#|      -acceptable values :
#|        -0 :-only the first seg will have N-tail and the last seg C-tail ;;
#|        -1 :-any segment (fragment) within each chain will have tails ;
#|      -default value :-0 ;;
  variable tails 2
  variable tailN $tails
  variable tailC $tails
  variable internalTails 0

#|    -'keepGlycine', 'keepGapGly', 'glyGapKeep' :
#|      -specifies whether glycine rec-gap residues should be ignored or
#|       _ mutated according the 'mutateGaps' arg .
#|      -acceptable values :-0 .-1 ;
#|      -default value :-1 ;;
#|    -'keepProline', 'keepGapPro', 'proGapKeep' :
#|      -specifies whether proline rec-gap residues should be ignored or
#|       _ mutated according the 'mutateGaps' arg .
#|      -acceptable values :-0 .-1 ;
#|      -default value :-1 ;;
#|    -'keepCysteine, 'keepGapCys', 'cysGapKepp' :
#|      -specifies whether cysteine rec-gap residues should be ignored or
#|       _ mutated according the 'mutateGaps' arg .
#|      -acceptable values :-0 .-1 ;
#|      -default value :-0 ;;
  variable keepGlycine 1
  variable keepProline 1
  variable keepCysteine 0

#|    -'resAtm', 'pivotTxt', 'pivotSelTxt', 'uniqueResAtom' :
#|      -selection text specifying one atom appearing once on each residue .
#|      -default value :-"name CA" ;;
#|    -'nTerPatch', 'firstPatch' :
#|      -N-terminal patch (PRES) name (psfgen keyword 'first') .
#|      -default value :-"NNEU" ;;
#|    -'cTerPatch' 'lastPatch' :
#|      -C-terminal patch (PRES) name (psfgen keyword 'last') .
#|      -default value :-"CNEU" ;;
#|    -'nTerGlyPatch', 'nTerPatchGly', 'firstPatchGly', 'firstGlyPatch' :
#|      -N-terminal patch (PRES) name for Gly (psfgen keyword 'first') .
#|      -default value :-"NGNE" ;;
#|    -'nTerProPatch', 'nTerPatchPro' 'firstPatchPro', 'firstProPatch' :
#|      -N-terminal patch (PRES) name for Pro (psfgen keyword 'first') .
#|      -default value :-"ACP" ;;
  variable resAtm "name CA"
  variable nTerPatch "NNEU"
  variable cTerPatch "CNEU"
  variable nTerGlyPatch "NGNE"
  variable nTerProPatch "ACP"

#|    -'chargeTarget', 'targetCharge', 'mutateChargeTarget' :
#|      -***incomplete implementation*** .
#|      -acceptable values :
#|        -"0", "zero", "neutral" :
#|          -a neutral charge for the whole system is sought .
#|          -two options will be implemented in future versions :
#|            -all residues neutral .
#|            -other strategies ;;
#|        -"" :-the chargeTarget arg will be ignored ;;
#|      -default value :-"" ;;
  variable chargeTarget ""

#|    -'bll', 'baseLogLevel', 'baseLogLvl' :
#|      -base log level used for log messages in the proc .
#|      -defines the values of ll1, ll2, and ll3 .
#|      -see library logLib .
#|      -defaul value :-1 ;;;
  variable bll 1

#|  -internal variables :- ;
  variable ll1 $bll
  variable ll2 [expr {$bll + 1}]
  variable ll3 [expr {$bll + 2}]
  variable singleMol 1
  variable id
  variable idL
  variable idR
  variable procN
  variable topCad ""
  variable resids
  variable pgnOut
  variable l_segR {}
  variable l_segL {}
  variable l_ligFAtInd {}
  variable selTxtLCad {}
  variable selTxtRCad {}
  variable comLine ""

#|  -procedures (namespace lr_trimComplex) :
#|    -proc ini_lr_trimComplex {} :
#|      -initializes internal namespace variables ;
  proc init_lr_trimComplex {} {
    global lr_trimComplex_version
    set_logName "lr_trimComplex"
    set_logVersion "$lr_trimComplex_version"
    }   ;# proc init_lr_trimComplex

#|    -proc setCLOptions {args} :
#|      -reads user options from command-line ;
  proc setCLOptions {args} {
# namespace variables
    variable src
    variable selId
    variable frame
    variable pdbId
    variable pdbIdL
    variable pdbIdR
    variable selTxtL
    variable selTxtR
    variable selTxtL_noInt
    variable l_topFile
    variable topDir
    variable l_topExt
    variable pdbAliasCad
    variable userPatchesCad
    variable prefix
    variable workPath
    variable pgnWrite
    variable segPrefix
    variable cutoff
    variable gapMax
    variable mutateGaps
    variable gapMut
    variable tail
    variable tailN
    variable tailC
    variable internalTails
    variable keepGlycine
    variable keepProline
    variable keepCysteine
    variable resAtm
    variable nTerPatch
    variable cTerPatch
    variable nTerGlyPatch
    variable nTerProPatch
    variable chargeTarget
    variable procN
    variable ll3
# internal variables
    set args_last {}
# read input command-line arguments (excluding logLib arguments)
    if {[expr {[llength $args]%2}] == 0} {
      logMsg " [get_logName_version]: processing list of arguments: $args" $ll3
      foreach {arg val} $args {
        logMsg " reading arg-val: $arg $val" $ll3
        switch [string tolower $arg] {
          "src" - "source" - "input" - "inputsource" {set src $val}
          "selid" - "selinfo" - "selinfoid" - "selinfokey" {set selId $val}
          "frame" - "frm" - "index" {set frame $val}
          "pdbid" - "pdbfile" - "id" - "vmdid" - "molid" {set pdbId $val}
          "pdbidl" - "pdbfilel" - "idl" - "vmdidl" - "molidl" {set pdbIdL $val}
          "pdbidr" - "pdbfiler" - "idr" - "vmdidr" - "molidr" {set pdbIdR $val}
          "seltxtl" - "l_seltxtl" - "seltxtls" {set selTxtL $val}
          "seltxtr" {set selTxtR $val}
          "seltxtl_noint" - "l_seltxtl_noint" - "seltxtl_noninteract" \
            - "seltxt_noninteractl" - "extralig" - "extraligands" \
            {set selTxtL_noInt $val}
          "l_topfile" - "l_topfiles" - "topfile" - "topfiles" \
            - "l_topologyfile" - "topologyfiles" {set l_topFile $val}
          "topdir" - "toppar" - "topologydir" - "topologypath" {set topDir $val}
          "l_topext" - "l_topfileext" {set l_topExt $val}
          "pdbaliascad" {set pdbAliasCad $val}
          "userpatchescad" - "userpatchcad" - "patchcad" - "patchescad" \
            - "patches" - "userpatches" {set userPatchesCad $val}
          "prefix" - "outprefix" - "outputprefix" {set prefix $val}
          "workpath" - "workfolder" - "workdir" - "outpath" - "outfolder" \
            - "outdir" {set workPath $val}
          "pgnwrite" - "pgnscript" - "pgnfile" - "psfgenfile" - "psfgenscript" \
            {set pgnWrite $val}
          "segprefix" - "segmentprefix" - "segletter" {set segPrefix $val}
          "pgnoutpath" - "pgndirout" - "pgndiroutput" {set pgnOutPath $val}
          "cutoff" - "distcutoff" - "distance" - "dist" - "r" {set cutoff $val}
          "gapmax" - "maxgap" - "resgap" - "gap" - "gaps" - "gapsize" - \
            "gaplen" - "lengap" - "gaplength" - "lengthgap" {set gapMax $val}
          "mutategaps" - "gapsmutate" - "mutate" {set mutateGaps $val}
          "gapmut" - "gapmutation" - "gapres" - "gapresidue" - "gapresname" \
            {set gapMut $val}
          "tails" - "tailslength" - "lengthtails" - "tailslen" - "lentails" \
            {set_tails $val}
          "tailn" - "ntail" - "tail-n" - "n-tail" {set tailN $val}
          "tailc" - "ctail" - "tail-c" - "c-tail" {set tailC $val}
          "internaltails" - "tailsinternal" - "intTails" \
            {set internalTails $val}
          "keepglycine" - "keepgapgly" - "glygapkeep" {set keepGlycine $val}
          "keepproline" - "keepgappro" - "progapkeep" {set keepProline $val}
          "keepcysteine" - "keepgapcys" - "cysgapkeep" {set keepCysteine $val}
          "resatm" - "pivottxt" -  "pivotseltxt" - "uniqueresatom" \
            {set resAtm $val}
          "nterpatch" - "firstpatch" {set nTerPatch $val}
          "cterpatch" - "lastpatch" {set cTerPatch $val}
          "nterglypatch" - "firstpatchgly" {set nTerGlyPatch $val}
          "nterpropatch" - "firstpatchpro" {set nTerProPatch $val}
          "chargetarget" - "targetcharge" - "mutatechargetarget" {}
          "bll" - "baseloglevel" - "baseloglvl" {set_bll $val}
          default {
            logMsg "$procN: argument (value) unknown: $arg ($val)" $bll
            set args_last [concat ${args_last} $arg $val]
            }
          }
        }
      logMsg "$procN: unused arguments: $args_last"
    } else {
      logMsg "$procN: Odd number of arguments: $args" $bll
      return ""
      }
    return $args_last
    }   ;# proc setCLOptions

#|    -proc setBaseLogLevel {} :
#|      -sets the values for variables ll1, ll2, and ll3 ;
  proc set_bll {val} {
    variable bll; variable ll1; variable ll2; variable ll3
    set bll $val; set ll1 $bll
    set ll2 [expr {$bll + 1}]; set ll3 [expr {$bll + 2}]
    logMsg " using base log level: $bll" $ll2
    logMsg " ll1, ll2, and ll3 set to: $ll1, $ll2, and $ll3, resp." $ll3
    }   ;# proc setBaseLogLevel

#|    -proc set_tails :
#|      -sets the values of tails, tailN, and tailC ;
  proc set_tails {val} {
    variable tails; variable tailN; variable tailC; variable ll3
    set tails $val
    set tailN $tails
    set tailC $tails
    logMsg " tails, tailN, and tailC values set to: $tails" $ll3
    }   ;# proc set_tails

#|    -proc loadMol {} :
#|      -load required molecules and configures them ;
  proc loadMol {complexType} {
    global selInfo
    variable ll1; variable ll2; variable ll3
    variable pdbId; variable pdbIdL; variable pdbIdR
    variable singleMol; variable idL; variable idR
    variable selId; variable selTxtL; variable selTxtR
    variable src; variable frame; variable workPath
    switch [string tolower $complexType] {
      "pdb" {   ;# single PDB file with both the lig and the rec
        set singleMol 1
        switch [string tolower $src] {
          "download" - "pdbid" - "pdbload" - "pdb" {
            set id [mol pdbload $pdbId]
            set idL $id; set idR $id
            logMsg " PDB model downloaded into id $id: $pdbId" $ll2
            animate goto $frame
            [atomselect $id "all" frame $frame] writepdb "${workPath}baseLR.tmp"
            logMsg "  written ${workPath}baseLR.tmp file" $ll3
            }
          "file" - "loadfile" - "localfile" - "localpdb" {
            set id [mol new $pdbId type pdb waitfor all]
            set idL $id; set idR $id
            logMsg " loaded local file into id $id: $pdbId" $ll2
            animate goto $frame
            set name [string range $pdbId \
              [expr {[string last / $pdbId] + 1}] \
              [expr {[string last ".pdb" $pdbId] - 1}]]
            mol rename $id $name
            logMsg " mol renamed to: $name" $ll3
            [atomselect $id "all" frame $frame] writepdb "b${workPath}aseLR.tmp"
            logMsg "  written ${workPath}baseLR.tmp file" $ll3
            }
          "id" - "ids" - "molid" - "molids" - "loaded" - "vmdid" - "vmd" {
            set id $pdbId
            set idL $id; set idR $id
            logMsg " using molecule ([molinfo $id get name]) VMD Id: $id" $ll2
            if {[string last ".pdb" [molinfo $id get name]] >= 0} {
              set name [string range [molinfo $id get name] 0 \
                [expr {[string last ".pdb" [molinfo $id get name]] - 1}]]
              mol rename $id $name
              logMsg " mol renamed to: $name" $ll3
              }
            [atomselect $id "all" frame $frame] writepdb "${workPath}baseLR.tmp"
            logMsg "  written ${workPath}baseLR.tmp file" $ll3
            }
          default {
            logMsg "Unknown option for 'src' arg: $src" $ll1
            return ""
            }
          }
        }
      "selinfo" {
        if {$selId != ""} {
          logMsg " reading input info from selInfo's selId: $selId" $ll2
          if {([info exists selInfo($selId,ligId)]) && \
              ([info exists selInfo($selId,recId)]) && \
              ([info exists selInfo($selId,ligSelTxt)]) && \
              ([info exists selInfo($selId,recSelTxt)]) && \
              ([info exists selInfo($selId,title)])} {
            # updating default selTxt* values
            if {$selTxtL == "not (protein or nucleic or water)"} {
              set selTxtLCad $selInfo($selId,ligSelTxt)
              }
            if {$selTxtR == "protein or nucleic"} {
              set selTxtRCad $selInfo($selId,recSelTxt)
              }
            logMsg " setting selTxtLCad: $selTxtLCad" $ll3
            logMsg " setting selTxtRCad: $selTxtRCad" $ll3
            if {[info exists selInfo($selId,frame)]} {
              set frame $selInfo($selId,frame)
              logMsg " setting frame: $frame" $ll3
              }
            if {$selInfo($selId,ligId) == $selInfo($selId,recId)} {
              set singleMol 1
              set id $selInfo($selId,ligId)
              set idL $id; set idR $id
              logMsg " setting id (same for both lig and rec): $id" $ll3
              [atomselect $id "all" frame $frame] writepdb "${workPath}baseLR.tmp"
              logMsg "  written ${workPath}baseLR.tmp file" $ll3
            } else {
              set singleMol 0
              set idL $selInfo($selId,ligId)
              set idR $selInfo($selId,recId)
              logMsg " setting idL: $idL" $ll3
              logMsg " setting idR: $idR" $ll3
              [atomselect $idL "all" frame $frame] writepdb "${workPath}baseL.tmp"
              logMsg "  written ${workPath}baseL.tmp file" $ll3
              [atomselect $idR "all" frame $frame] writepdb "${workPath}baseR.tmp"
              logMsg "  written ${workPath}baseR.tmp file" $ll3
              }
          } else {
            logMsg "$procN: selInfo missing data for selId: $selId"
            logMsg ", required array keys selInfo($selId,<key>): ligId, recId, ligSelTxt, recSelTxt, and title" $ll1
            return ""
            }
        } else {
          logMsg "the 'selInfo' option requires the 'selId' arg." $ll1
          return ""
          }
        }
      "dlg" {   ;# lig and rec from different mols including a .dlg file
        # not implemented yet
        }
      "-h" - "--help" {   ;# print help
        lr_trimComplex_help
        return ""
        }
      default {
        logMsg "Unknown option for 'complexType' arg: $src" $ll1
        return ""
        }
      }
    }   ;# proc loadMol

#|    -proc atomSel {} :
#|      -create atom selections ;
    proc atomSel {} {
      variable ll1; variable ll2; variable ll3
      variable singleMol; variable selTxtL_noInt
      variable frame; variable resAtm; variable selTxtLCad
      variable cutoff; variable selTxtL; variable selTxtR; variable selTxtRCad
      variable idL; variable idR; variable l_segL; variable l_ligFAtInd
      if {$singleMol} {
        set ligSel ""; set interactCASel ""   ;# provisional
        if {($ligSel == "") || ($interactCASel == "")} {
          logMsg " using interaction cutoff value: $cutoff" $ll2
          logMsg " using selTxtLs: $selTxtL" $ll2
          logMsg " using selTxtR: $selTxtR" $ll2
          logMsg " using selTxtL_noInt: $selTxtL_noInt" $ll2
          logMsg " using resAtm: $resAtm" $ll2
          set l_ligFAtInd {}
          set l_segL {}
          set i 1
#          set l_selTxtL_prov $selTxtL
          set selTxtLCad ""
          foreach sTxt ${selTxtL} {
            set ligSel [atomselect $idR "$sTxt"]   ;# frame needed?
            if {$selTxtLCad == ""} {set selTxtLCad "($sTxt)"
              } else {set selTxtLCad "$selTxtLCad or ($sTxt)" }
            lappend l_segL "L${i}"
            lappend l_ligFAtInd [lindex [$ligSel get index] 0]
            incr i
            $ligSel delete
            }
          logMsg " list of rec-interacting ligands segNames: $l_segL" $ll3
          logMsg " selTxtLCad: $selTxtLCad"
          foreach sTxt ${selTxtL_noInt} {
            set ligSel [atomselect $idR "$sTxt"]   ;# note: no frm considered
            lappend l_segL "L${i}"
            lappend l_ligFAtInd [lindex [$ligSel get index] 0]
            incr i
            $ligSel delete
            }
          if {${selTxtL_noInt} != ""} {
            logMsg " amended list of ligands segNames: $l_segL" $ll3}
          set selTxtRCad "(same residue as protein within $cutoff of ($selTxtLCad)) and ($selTxtR) and ($resAtm)"
          logMsg " using interactCASel selTxt: $selTxtRCad" $ll2
          set interactCASel [atomselect $idR "$selTxtRCad"]
        } else {
          logMsg " atom selections previously defined:" $ll2
          logMsg "  lig: [$ligSel text] " $ll2
          logMsg "  rec: [$interactCASel text]" $ll2  
          }
        logMsg " rec residues selected: \
                [$interactCASel get {chain resname resid}]" $ll2
      } else {
        # not implemented yet
        logMsg " lig and rec from different molecules not yet implemeted" $ll1
        return ""
        }
#      $ligSel delete
      $interactCASel delete
      }   ;# proc atomSel

#|    -proc outputConfig {} :
#|      -check workPath, prefix, and pgnWrite arguments ;
    proc outputConfig {} {
      variable ll1; variable ll2; variable ll3
      variable workPath; variable prefix; variable idL; variable idR
      variable pgnOut; variable pgnWrite; variable selTxtL; variable selTxtR
      variable frame; variable cutoff; variable selTxtLCad; variable selTxtRCad
      variable pgnOutPath
# 
      set ligSel [atomselect $idL $selTxtLCad]
      set interactCASel [atomselect $idR $selTxtRCad]
      switch [string tolower $workPath] {
        "" - "." - "./" {
          set workPath "./"
          logMsg "using current directory as working path: [pwd]" $ll2
          }
        default {
          if {[string index $workPath end] != "/"} {
            set workPath "$workPath/"
            logMsg "character appended to the workPath '/'" $ll3
            }
          if {![file exist $workPath]} {
            exec mkdir -p $workPath
            logMsg "working path no found; now created: $workpath" $ll2
            }
          }
        }
      if {($prefix == "") || ($prefix == "auto")} {
        set prefix "[molinfo $idR get name]-[lindex [$ligSel get resname] 0]_${cutoff}"
        logMsg " generated prefix: $prefix" $ll3
        }
      logMsg " setting prefix for output files: $prefix" $ll2
      if {$pgnWrite} {
        set pgnOut [open "${workPath}${prefix}.pgn.tcl" w]
        logMsg " pgn script for psfgen created: ${workPath}${prefix}.pgn.tcl" $ll2
      } else {
        logMsg " pgn script for psfgen will not be written" $ll2
        }
      # create dir for the psfgen output if it does not exists
      if {($pgnOutPath != "") || ($pgnOutPath != ".") || ($pgnOutPath != "./")} {
        logMsg "using pgnOutPath: $pgnOutPath" $ll3
        if {![file exists $pgnOutPath]} {
          exec mkdir -p $pgnOutPath
          logMsg "  created dir for psfgen output: $pgnOutPath" $ll3
          }
        }
      $ligSel delete
      $interactCASel delete
      }   ;# proc outputConfig

#|    -proc readTopology {} :
#|      -read topology files ;
    proc readTopology {} {
      variable ll1; variable ll2; variable ll3
      variable l_topFile; variable topDir; variable topCad
      if {($l_topFile == {}) && ($topDir == "")} {
        logMsg "'l_topFile' and/or 'topDir' arguments must be specified" $ll1
        return ""
        }
      set topCad ""
      if {[llength $l_topFile] >= 1} {
        logMsg " reading topology files from l_topFile: $l_topFile" $ll2
        foreach topFile $l_topFile {
          if {[file exists $topFile]} {
            logMsg " topology file found: $topFile" $ll3
            set topCad "${topCad}topology $topFile\n"
          } else {
            logMsg " topology file not found: $topFile" $ll2
            }
          }
        }
      if {$topDir != ""} {
        foreach dir $topDir {
          logMsg " searching top files (ext: $l_topExt) in dir: $dir" $ll2
          foreach ext $l_topExt {
            if {[catch {lsort [eval glob "${dir}*.${ext}"]} topDirFiles]} {
              set topDirFiles ""
              logMsg " no .${ext} files found in dir ${dir}" $ll3
              }
            foreach topFile $topDirFiles {
              logMsg " topology file found: $topFile" $ll3
              set topCad "${topCad}topology $topFile\n"
              }
            }
          }
        }
      if {$topCad == ""} {
        logMsg "no topology files were found" $ll1
        return ""
        }
      }   ;# proc readTop

#|    -proc detectChains :
#|      -unused in version 004 .
#|      -pending to be fixed ;
    proc detectChains {} {
      variable ll1; variable ll2; variable ll3
      variable singleMol; variable resids
      variable selTxtR; variable selTxtRCad; variable idR
 
      if {$singleMol} {

     } else {
        infoMsg "chain detection not implemented yet for multiple mols" $ll1
        return ""
        }
      }   ;# proc detectChains

#|    -proc pgnWriteHeader {} :
#|      -header for pgn file ;
    proc pgnWriteHeader {} {
      variable ll1; variable ll2; variable ll3; variable comLine
      variable topCad; variable pdbAliasCad; variable pgnWrite; variable pgnOut
      if {($topCad == "") || ($pdbAliasCad == "")} {
        logMsg "both topCad and pdbAliadCad must be defined" $ll1
        return ""
        } else {
          if {$pgnWrite} {
            puts $pgnOut "# psfgen pgn file created by: [get_logName_version] on [clock format [clock seconds]]"
            puts $pgnOut "# vmd command: $comLine"
            puts $pgnOut ""
            puts $pgnOut $topCad
            puts $pgnOut $pdbAliasCad
          } else {
#            package require psfgen
#            eval $topCad
#            eval $pdbAliasCad
            }
          logMsg " psfgen header: topCad:" $ll2
          logMsg $topCad $ll2
          logMsg " psfgen header: pdbAliasCad:" $ll2
          logMsg $pdbAliasCad $ll2
          }
      }   ;# proc pgnWriteHeader

#|    -proc addRes_firstTail :
#|      -defines a new segment .
#|      -adds a residue and first patch .
#|      -tries to mutate it ;
    proc addRes_firstTail {resId resName chain segName tailRes} {
      variable ll1; variable ll2; variable ll3
      variable pgnWrite; variable pgnOut
      variable mutateGaps; variable gapMut;
      variable nTerPatch; variable nTerGlyPatch; variable nTerProPatch
      variable keepGlycine; variable keepProline; variable keepCysteine

      set i $resId
      set seg $segName
      if {$pgnWrite} {
        puts $pgnOut "segment ${seg} \{"
        puts $pgnOut "  residue $i $resName $chain"
      } else {
#        eval "segment ${seg} \{"
#        eval "  residue $i $resName $chain"
        }
      # mutate residue?
      switch $resName {
        "GLY" {
          if {$mutateGaps && !$keepGlycine && $tailRes} {
            logMsg " mutating first residue in segment $seg to $gapMut"
            logMsg " first patch for segment $seg: $nTerPatch" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  first $nTerPatch"
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "first $nTerPatch"
#              eval "mutate $i $gapMut"
              }
          } else {
            logMsg " first patch for segment ${seg}: $nTerGlyPatch" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  first $nTerGlyPatch"
            } else {
#              eval "first $nTerGlyPatch"
              }
            }
          }
        "PRO" {
          if {$mutateGaps && !$keepProline && $tailRes} {
            logMsg " mutating first residue in segment $seg to $gapMut"
            logMsg " first patch for segment $seg: $nTerPatch" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  first $nTerPatch"
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "first $nTerPatch"
#              eval "mutate $i $gapMut"
              }
          } else {
            logMsg " first patch for segment ${seg}: $nTerProPatch" $ll3 
            if {$pgnWrite} {
              puts $pgnOut "  first $nTerProPatch"
            } else {
#              eval "first $nTerProPatch"
              }
            }
          }
        "CYS" {
          if {$mutateGaps && !$keepCysteine && $tailRes} {
            logMsg " mutating first residue in segment $seg to $gapMut"
            logMsg " first patch for segment $seg: $nTerPatch" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  first $nTerPatch"
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "first $nTerPatch"
#              eval "mutate $i $gapMut"
              }
          } else {
            logMsg " first patch for segment $seg: $nTerPatch" $ll3 
            if {$pgnWrite} {
              puts $pgnOut "  first $nTerPatch"
            } else {
#              eval "first $nTerPatch"
              }
            }
          }
        "ALA" {
          logMsg " first patch for segment $seg: $nTerPatch" $ll3
          if {$pgnWrite} {
            puts $pgnOut "  first $nTerPatch"
          } else {
#            eval "first $nTerPatch"
            }
          }
        default {
          if {$mutateGaps && $tailRes} {
            logMsg " mutating first residue in segment $seg to $gapMut"
            if {$pgnWrite} {
              puts $pgnOut "  first $nTerPatch"
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "first $nTerPatch"
#              eval "mutate $i $gapMut"
              }
          } else {
            if {$pgnWrite} {
              puts $pgnOut "  first $nTerPatch"
            } else {
#              eval "first $nTerPatch"
              }
            }
          logMsg " first patch for segment $seg: $nTerPatch" $ll3
          }
        }


      }   ;# proc addRes_firstTail


#|    -proc addRes_gap :
#|      -add a residue and tries to mutate it ;
    proc addRes_gap {resId resName chain segName} {
      variable ll1; variable ll2; variable ll3
      variable pgnWrite; variable pgnOut
      variable mutateGaps; variable gapMut;
      variable keepGlycine; variable keepProline; variable keepCysteine
      set i $resId
      set seg $segName
      logMsg " adding residue $resName $i to segment $seg" $ll3
      if {$pgnWrite} {
        puts $pgnOut "  residue $i $resName $chain"
      } else {
#        eval "residue $i $resName $chain"
        }
      # mutate gap resid?
      switch $resName {
        "GLY" {
          if {$mutateGaps && !$keepGlycine} {
            logMsg "GLY $i mutated to $gapMut in seg $segName" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "  mutate $i $gapMut"
              }
            }
          }
        "PRO" {
          if {$mutateGaps && !$keepProline} {
            logMsg "PRO $i mutated to $gapMut in seg $segName" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "  mutate $i $gapMut"
              }
            }
          }
        "CYS" {
          if {$mutateGaps && !$keepCysteine} {
            logMsg "CYS $i mutated to $gapMut in seg $segName" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "  mutate $i $gapMut"
              }
            }
          }
        "ALA" {logMsg "gap ALA $i left intact" $ll3}
        default {
          if {$mutateGaps} {
            logMsg "$resName $i mutated to $gapMut in seg $seg" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "  mutate $i $gapMut"
              }
            }
          }
        }
      }   ;# proc addRes_gap


#|    -proc addRes_lastTail {resId resName chain segName tailRes} :
#|      -ends a segment .
#|      -adds a residue and last patch .
#|      -tries to mutate it ;
    proc addRes_lastTail {resId resName chain segName tailRes} {
      variable ll1; variable ll2; variable ll3
      variable pgnWrite; variable pgnOut
      variable mutateGaps; variable gapMut;
      variable cTerPach
      variable keepGlycine; variable keepProline; variable keepCysteine

      set i $resId
      logMsg " added resid $resName $i as C-tail" $ll3
      if {$pgnWrite} {
        puts $pgnOut "  residue $i $resName $chain"
      } else {
#       eval "  residue $i $resName $chain"
        }
      # mutate tail residue?
      switch $resName {
        "GLY" {
          if {$mutateGaps && !$keepGlycine && $tailRes} {
            logMsg " C-tail residue GLY $i mutated to $gapMut" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "  mutate $i $gapMut"}
            }
          }
        "PRO" {
          if {$mutateGaps && !$keepProline && $tailRes} {
            logMsg " C-tail residue PRO $i mutated to $gapMut" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "  mutate $i $gapMut"
              }
            }
          }
        "CYS" {
          if {$mutateGaps && !$keepCysteine && $tailRes} {
            logMsg " C-tail residue CYS $i mutated to $gapMut" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "  mutate $i $gapMut"
              }
            }
          }
        default {
          if {$mutateGaps && ($resName != $gapMut)} {
            logMsg " C-tail residue $resName $i mutated to $gapMut" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  mutate $i $gapMut"
            } else {
#              eval "  mutate $i $gapMut"
              }
            }
          }
        }
      }   ;# proc addRes_lastTail

#|    -proc addRes_lig {} :
#|      -adding segments for ligand residues ;
    proc addRes_lig {} {
      variable ll1; variable ll2; variable ll3
      variable l_segR; variable pgnWrite; variable pgnOut; variable singleMol
      variable workPath; variable prefix; variable l_segL; variable selTxtL
      variable idL; variable l_ligFAtInd; variable selTxtLCad

      logMsg " adding ligand segments" $ll3
      if {$singleMol} {set suf "LR"} else {set suf "R"}
      if {$pgnWrite} {
        set i 0
        foreach segL $l_segL {
          set ligSel [atomselect $idL "index [lindex $l_ligFAtInd $i]"]
          set resLCad [$ligSel get {resid resname chain}]
          logMsg "ligand residue added to seg $segL: [lindex $resLCad 0]" $ll2
          puts $pgnOut "segment $segL \{"
          puts $pgnOut "  first NONE"
          puts $pgnOut "  last NONE"
          puts $pgnOut "  residue [lindex $resLCad 0]"
          puts $pgnOut "  \}\n"
          incr i
          $ligSel delete
          }
      } else {
        set i 0
        foreach segL $l_segL {
          set ligSel [atomselect $idL "index [lindex $l_ligFAtInd $i]"]
          set resLCad [$ligSel get {resid resname chain}]
          logMsg "ligand residue added to seg $segL: $resLCad" $ll2
#          eval "segment $segL \{"
#          eval "  first NONE"
#          eval "  last NONE"
#          eval "  residue $resLCad"
#          eval "  \}"
          incr i
          $ligSel delete
          }
        }
      }   ;# proc addRes_lig


#|    -proc pgnWriteTail {} :
#|      -adding ligands and final keywords in the pgn file ;
    proc pgnWriteTail {} {
      variable ll1; variable ll2; variable ll3
      variable l_segL; variable l_segR; variable pgnWrite; variable pgnOut
      variable singleMol; variable workPath; variable prefix
      variable userPatchesCad; variable pgnOutPath

      if {$singleMol} {set suf "LR"} else {set suf "R"}
      if {$pgnWrite} {
        # add coordpdb keywords
        foreach seg $l_segR {
          puts $pgnOut "coordpdb base${suf}.tmp $seg"
          }
        foreach seg $l_segL {
          puts $pgnOut "coordpdb base${suf}.tmp $seg"
          }
        # add last keywords
        if {$userPatchesCad == ""} {
          puts $pgnOut ""
          logMsg " user patch(es) string not included" $ll3
        } else {
          puts $pgnOut "\n$userPatchesCad\n"
          logMsg " user patch(es) added to pgnFile:\n$userPatchesCad" $ll2
          }
        puts $pgnOut "guesscoord\n"
        puts $pgnOut "writepsf ${pgnOutPath}${prefix}.psf"
        puts $pgnOut "writepdb ${pgnOutPath}${prefix}.pdb\n"
      } else {
        # add coordpdb keywords
        foreach seg $l_segR {
#          eval "coordpdb base${suf}.tmp $seg"
          }
        foreach seg $l_segL {
#          eval "coordpdb base${suf}.tmp $seg"
          }
        # add last keywords
#        eval "$regenerateCad\n"
#        eval "guesscoord\n"
#        eval "writepsf ${workPath}${prefix}.psf"
#        eval "writepdb ${workPath}${prefix}.pdb\n"
        }
      }   ;# proc pgnWriteTail

#|    - ;

#|  -proc lr_trimComplex { } :
#|    -writes pdb files for each fragment in the protein that interacts with
#|     _ the ligand .
#|    -notes :
#|      -topology file/dir specification mandatory (no default values) ;;

proc lr_trimComplex {complexType args} {
# global variables
  global selInfo

# namespace variables
  variable ll1; variable ll2; variable ll3
  variable procN; variable src; variable selId; variable frame; variable pdbId
  variable pdbIdL; variable pdbIdR; variable selTxtL
  variable selTxtR; variable selTxtRCad
  variable l_topFile; variable topDir; variable comLine; variable selTxtLCad
  variable l_topExt; variable pdbAliasCad; variable prefix; variable workPath
  variable pgnWrite; variable segPrefix; variable cutoff; variable gapMax
  variable mutateGaps; variable gapMut; variable tail; variable tailN
  variable tailC; variable internalTails; variable keepGlycine
  variable keepProline; variable keepCysteine; variable resAtm
  variable nTerPatch; variable cTerPatch; variable nTerGlyPatch
  variable nTerProPatch; variable chargeTarget; variable pgnOut
  variable singleMol; variable idL; variable idR; variable l_segR
 
# proc name
  set procN [lindex [info level 0] 0]

# interpreting logLib and namespace arguments 
  set args_rest [eval arg_interpreter $args]
  logMsg " arg_rest after logLib arg_interp: $args_rest" $ll3
  set args_rest [eval setCLOptions $args_rest]

# initial info messages
  logMsg "+++ Trim of ligand-receptor complexes +++" $ll1
  set comLine "$procN $complexType $args"
  logMsg " command line: $comLine" $ll1
  logMsg " unused arguments: $args_rest" $ll1

# loadMolecule(s)
  loadMol $complexType

# create atom selections
  atomSel

# configure output
  outputConfig

# read toplogy files
  readTopology

# detect chain fragments

        set interactCASel [atomselect $idR $selTxtRCad]
        foreach chain [$interactCASel get chain] rid [$interactCASel get resid] {
          if {[info exists resids($chain)]} {
            lappend resids($chain) $rid
          } else {
            set resids($chain) $rid
            }
          }
        logMsg " chains in contact with the ligand(s): [array names resids]" $ll2
        $interactCASel delete
        set interactCASel ""
 

  logMsg " using parameters for fragment building:" $ll2
  logMsg "  gapMax: $gapMax  tailN: $tailN  tailC: $tailC" $ll2
  logMsg "  mutateGaps: $mutateGaps  gapMutation: $gapMut" $ll2
  logMsg "  keepProline: $keepProline  keepCysteine: $keepCysteine" $ll2
  logMsg "  chargeTarget: $chargeTarget" $ll3
  logMsg "  NTerPatch: $nTerPatch  CTerPatch: $cTerPatch" $ll2
  logMsg "  NTerGlyPatch: $nTerGlyPatch  NTerProPatch: $nTerProPatch" $ll2
  if {$tailN > [expr {$gapMax + 1}]} {
    set tailN [expr {$gapMax + 1}]
    logMsg " 'tailN' value adjusted to 'gapMax' + 1: $tailN" $ll2
    }
  if {$tailC > [expr {$gapMax + 1}]} {
    set tailC [expr {$gapMax + 1}]
    logMsg " 'tailC' value adjusted to 'gapMax' + 1: $tailN" $ll2
    }

# header for pgn file
  pgnWriteHeader

  if {$singleMol} {
    # process each chain
    set ll_frag_rid {}
    set l_segR {}
    foreach chain [array names resids] {
      logMsg " processing chain: $chain" $ll2
      set prevrid 0   ;# last resid assigned to a fragment (within a chain)
      set l_rid {}   ;# list of resids for the current fragment
      set l_gapTail {}   ;# list of resids added to the current gap
      set nFrag 1
      set r 0
      while {$r < [llength $resids($chain)]} {
        set rid [lindex $resids($chain) $r]
        set seg "${chain}${nFrag}"
        logMsg " processing contact resid for segment $seg: $rid" $ll3
        # filling a gap or adding an N-tail
        if {[llength $l_rid] == 0} {   ;# first interacting resid
          # adds an N-tail if there are residues available
          set i [expr {$rid-$tailN}]
          while {$i <= $rid} {   ;#
            if {!$internalTails && ($nFrag > 1)} {   ;# no N-tail in segment
              set i $rid   ;# makes this loop run only once
              logMsg " avoided adding N-tail to segment: $seg" $ll3
              }
            set tailId [atomselect $idR "chain $chain and resid $i and $resAtm"]
            if {[$tailId num] == 1} {
              # start psfgen segment and add first residue
              logMsg " starting segment ${seg}" $ll3
              logMsg " added resid [$tailId get resname] $i as N-tail" $ll3
              addRes_firstTail $i [$tailId get resname] $chain $seg [expr {$i != $rid}]
              if {$i == $rid} {incr r}
              lappend l_rid $i
              lappend l_segR $seg
              $tailId delete
              break
            } else {
              logMsg "unable to add resid $i as N-tail" $ll2
              }
            $tailId delete
            incr i
            }   ;# while
        } else {   ;# next interacting resid within a chain/fragment
# *** check adding first interaction sphere resids***
          set prevrid [lindex $l_rid end]
          if {$rid == [expr {$prevrid+1}]} {  ;# consecutive resid
            set seqId [atomselect $idR "chain $chain and resid $rid and $resAtm"]
            logMsg " adding residue to segment ${seg}: [$seqId get resname] $rid" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  residue $rid [$seqId get resname] $chain"
            } else {
#              eval "residue $rid [$seqId get resname] $chain"
              }
            $seqId delete
            lappend l_rid $rid
            incr r
          } else {   ;# non-consecutive resid
            if {$prevrid >= [expr {$rid-$gapMax-1}]} {   ;# add gap resids
              for {set i [expr {$prevrid+1}]} {$i < $rid} {incr i} {
                set gapId [atomselect $idR "chain $chain and resid $i and $resAtm"]
                if {[$gapId num] == 1} {
                  addRes_gap $i [$gapId get resname] $chain $seg
                  lappend l_rid $i
                } else {
                  logMsg "missing residue in chain $chain: $i" $ll1
                  return ""
                  }
                $gapId delete
                }
            } else {   ;# new segment needed
              # adding segment C-tail
              logMsg " adding a new segment..." $ll2
              if {$internalTails} {
                logMsg "  adding internal C-tail for segment $seg"
                for {set i [expr {$prevrid+1}]} {$i <= [expr {$prevrid+$tailC}]} {incr i} {
                  set tailId [atomselect $idR "chain $chain and resid $i and $resAtm"]
                  if {[$tailId num] == 1} {
                    addRes_lastTail $i [$tailId get resname] $chain $seg [expr {$i != $rid}]
                    lappend l_rid $i
                  } else {
                    logMsg " unable to add residue $i tail of segnment; $seg" $ll1
                    break
                    }
                  $tailId delete
                  }
                }
              # adding 'last' patch and ending segment
              if {$pgnWrite} {
                puts $pgnOut "  last $cTerPatch"
                puts $pgnOut "  \}\n"
              } else {
#                eval "last $cTerPatch"
#                eval "\}\n"
                }
              logMsg " resids for complete segment ${seg}: $l_rid" $ll2
              lappend ll_frag_rid $l_rid
              logMsg " current ll_frag_rid: ${ll_frag_rid}" $ll2
              set l_rid {}
              incr nFrag
              }
            }
          }   ;# if (next interacting) 
        }   ;# while
        logMsg "  adding chain C-tail for segment $seg"
        set prevrid [lindex $l_rid end]
        for {set i [expr {$prevrid+1}]} {$i <= [expr {$prevrid+$tailC}]} {incr i} {
          set tailId [atomselect $idR "chain $chain and resid $i and $resAtm"]
          if {[$tailId num] == 1} {
            addRes_lastTail $i [$tailId get resname] $chain $seg [expr {$i != $rid}]
            lappend l_rid $i
          } else {
            logMsg " unable to add residue $i tail of segnment; $seg" $ll1
            break
            }
          $tailId delete
          }
          # adding 'last' patch and ending segment
          if {$pgnWrite} {
            puts $pgnOut "  last $cTerPatch"
            puts $pgnOut "  \}\n"
          } else {
#            eval "last $cTerPatch"
#            eval "\}\n"
            }
          logMsg " resids for complete segment ${seg}: $l_rid" $ll2
          lappend ll_frag_rid $l_rid
          logMsg " current ll_frag_rid: ${ll_frag_rid}" $ll2
      }   ;# foreach chain
  } else {
    logMsg "chain fragment detection not yet implemented for two molecs " $ll1
    return ""
    }
  # add ligand segment(s)
  addRes_lig
  # end pgn file
  pgnWriteTail

  if {$pgnWrite} {close $pgnOut}

  }   ;# proc lr_trimComplex

#|   -proc apply_restrictions {} :
#|     -will write pdf-format files specifying restrictions for namd2 calcs .
#|     -to be used after psfgen was used to run the script created by
#|      _ lr_trimComplex .
#|     -arguments :
#|       -psf :
#|         -.psf file created by psfgen .
#|         -acceptable values :
#|           -path and files name of a .psf file of a trimmed complex .
#|           -"", "auto" :
#|             -the name used in the lr_trimComplex run is used ;;
#|         -default value :-"auto" ;;
#|       -pdb :
#|         -.pdb file created by psfgen .
#|         -acceptable values :
#|           -path and files name of a .pdb file of a trimmed complex .
#|           -"", "auto" :
#|             -the name used in the lr_trimComplex run is used ;;
#|         -default value :-"auto" ;;
#|       -l_restrictType :
#|         -list of restrictions to be applied .
#|         -acceptable values :
#|           -"equil" :-***not implemented yet*** ;
#|           -"QM/MM", "QM-MM", "QMMM" :
#|             -interacting ligands will be unrestricted .
#|             -interacting residues will unrestricted .
#|             -non-interacting ligands will have harmonic restrictions .
#|             -noninteracting residue's sidechains will have harmonic
#|              _ armonic restrictions .
#|             -all backbone atoms will have harmonic restrictions .
#|             -the QM region will be the interacing ligs and residues ;;
#|         -default value :-"all" ;;;;
  proc apply_restrictions {{psfFile "auto"} {pdbFile "auto"} \
                           {l_restrictType "QM/MM"}} {
    # namespace variables:
    variable ll1; variable ll2; variable ll3
    variable selTxtR; variable selTxtL; variable selTxtL_noInt
    variable l_ligFAtInd; variable selTxtLCad; variable selTxtRCad
    variable workPath; variable prefix; variable pgnOutPath
    # load files
    foreach ext {"psf" "pdb"} {
      if {([set ${ext}File] == "auto") || ([set ${ext}File] == "")} {
        set ${ext}File "${pgnOutPath}${prefix}.${ext}"
        if {![file exists [set ${ext}File]]} {
          logMsg "File [set ${ext}File] not found." $ll1
          return
          }
        }
      }
    logMsg "\nlr_trim_complex::apply_restrictions:" $ll1
    package require topotools
    logMsg " list of restrictions: $l_restrictType" $ll1
    logMsg " psf file: $psfFile"
    logMsg " pdb file: $pdbFile"
    set idC [mol new $psfFile]
    mol addfile $pdbFile
    set syst [atomselect $idC "all"]
    logMsg " complex loaded: [$syst num] atoms" $ll2
    logMsg " total charge: [format "%.3f" [eval vecadd [$syst get charge]]]"
    # process types of restictions to be applied
    foreach restrictType $l_restrictType {
      switch [string tolower $restrictType] {
        "qm/mm" {
          # processing list of non-interacting ligands
          set noIntLigCad ""
          if {${selTxtL_noInt} != ""} {
            foreach sTxt ${selTxtL_noInt} {
              if {$noIntLigCad == ""} {set noIntLigCad "($sTxt)"
                } else {set noIntLigCad "$noIntLigCad or ($sTxt)"}
              }
            logMsg "noIntLigCad created: $noIntLigCad" $ll3
          } else {set noIntLigCad "none"}
          # applying restrictions
          set interactCASel [atomselect $idC "$selTxtRCad"]
          set harmRSel [atomselect $idC "(backbone or name HA HA1 HA2 HN) or (protein and not (resid [$interactCASel get resid] and sidechain)) or ($noIntLigCad)"]
          logMsg " selTxt for atoms with harmonic-potential restrictions: [$harmRSel text]" $ll2
          $syst set beta 0
          $syst set occupancy 0
          $harmRSel set beta 1
          $syst writepdb "${pgnOutPath}${prefix}_noInt.Bflag"
          logMsg " Bflag file for harmonic restrictions written: ${pgnOutPath}${prefix}_noInt.Bflag " $ll1
          # selecting QM region
          $interactCASel delete
          set interactCASel [atomselect $idC "($selTxtRCad) and not (resname GLY PRO)"]
          set qmReg [atomselect $idC "(protein and resid [$interactCASel get resid] and sidechain) or ($selTxtLCad)"]
          logMsg " selection text used for QM region: [$qmReg text]" $ll2
          set mmReg [atomselect $idC "(protein and resid [$interactCASel get resid] and name CA CB) or (protein and resid [$interactCASel get resid] and resname PRO and name N CD)"]
          topo guessatom element mass
          $syst set beta 0
          $syst set occupancy 0
          $qmReg set beta 1
          $mmReg set occupancy 1
          $syst writepdb "${pgnOutPath}${prefix}_qmreg.BOflag"
          logMsg " BOflag file for QM region written: ${pgnOutPath}${prefix}_qmreg.BOflag" $ll1
          logMsg " the QM region contains [$qmReg num] atoms" $ll2
          logMsg " the MM region contains [$mmReg num] atoms" $ll2
          # deleting atom selections
          $interactCASel delete
          $harmRSel delete
          $qmReg delete
          $mmReg delete
          }
        "all" {}
        default {}
        }
      }
    $syst delete
    }   ;# proc apply_restrictions

#|  -proc lr_trimComplex_help {{opt ""}} :- ;
proc lr_trimComplex_help {{opt ""}} {
  puts "+++ Trim of ligand-receptor complexes +++"
  puts "  usage: lr_trimComplex <complexType> \[<opt-arg value> ...\]"
  puts "  help option: $opt"
  switch [string tolower $opt] {
    "complextype" {
      puts "  complexType options: pdb"
      }
    "loglib" - "loglib.tcl" {
      ::logLib::logLib_help $opt
      }
    default {
      puts "  complexType options: pdb"
      puts "  optional args: src "
      puts "                 pdbId pdbIdL pdbIdR"
      puts "                 workPath prefix"
      puts "                 selTxtL selTxtR"
      puts "                 atmSelL atmSelR"
      puts "                 cutoff"
      puts "                 gapMax tails tailN tailC"
      puts "                 mutateGaps keepProline"
      puts "                 bll"
      puts "  required libraries: logLib.tcl regVar.tcl"
      puts "  logLib::logLib_help:"
      ::logLib::logLib_help $opt
      }
    }
  }   ;# proc lr_trimComplex_help 


#|  -proc readTop {} :- ;
proc readTop {} {
  }   ;# proc readTop

proc set_contactSelCAs {} {
  }   ;# proc set_contactSelCAs {}

#|  -proc genFrags {res } :
#|    -arguments :
#|      -res :
#|        -atomselectios with receptor alpha carbons for residues in
#|         _ direct contact with the ligand .
#|        -the selection may include several chains ;;;
proc genChainFrags {res} {

# classify resid by chain
  foreach chain [$res get chain] rid [$res get resid] {
    if {[info exists resids($chain)]} {
      lappend resids($chain) $rid
    } else {
      set resids($chain) $rid
      }
    }
  foreach chain [array names resids] {
    foreach rid $resids($chain) {
      puts "chain: $chain   resid: $rid"
      }
    }


# process fragments
  set nResFinal 1
  set l_gaps {}
  set l_rid {}
  set l_fragGaps {}

  }   ;# proc genChainFragments

  }   ;# namespace eval lr_trimComplex

::lr_trimComplex::init_lr_trimComplex

#|  - ;
#|- ;





