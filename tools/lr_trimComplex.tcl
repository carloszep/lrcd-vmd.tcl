#|-lr_trimComplex.tcl :
#|  -generate simplified versions of ligand-protein complexes .
#|  -generates inputs for psfgen and namd .
#|  -based on the old script extAChEFrag.tcl .
#|  -author :-Carlos Z. Gómez Castro ;
#|  -date :
#|    -created :-2025-01-13.Mon ;
#|    -modified :-2025-01-23.Thu ;;
#|  -version :
#|    -002 :
#|      -implementing new code for variable gap lengths ;
#|    -001 :
#|      -procedure's command-line interface .
#|      -procedure first defined .
#|      -original code tested ;
#|    -to do :
#|      -manage crystallographic water .
#|      -move to a namespace scheme ;;
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
#|          -"selInfo" :
#|            -the selection text and ids for lig and rec molecules are taken
#|             _ form the global selInfo array .
#|            -the 'selId' arg must be specified .
#|            -the complex may come from a single molecule or from two ;
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
#|            -"download", "pdbId", "pdbLoad", "PDB" :
#|              -pdb file(s) downloaded using the VMD command 'pdbload' .
#|              -arg 'pdbid' will be interpreted as a 4-letter PDBID ;
#|            -"file", "loadFile", "localFile", "localPDB" :
#|              -file names to locally load the files are assumed .
#|              -can contain absolute or relative file paths ;
#|            -"id", "ids", "molId", "molIds", "loaded", "vmdId", "vmd" :
#|              -the molecules are already loaded into VMD .
#|              -the VMD id(s) should be specified in args 'pdbid', 'pdbidl',
#|               _ or 'pdbidr' ;;
#|          -default value :-"id" ;;
#|        -'selId', 'selInfo', 'selInfoId', 'selInfoKey' :
#|          -selId key to use the information stored in the selInfo array .
#|          -implies that the lig and rec molecules were
#|           _ already loaded into VMD (maybe through a molInfo file) .
#|          -the atom selections registered in the selInfo array are expected
#|           _ to be more general (i.e., "protein and chain A", "not protein") .
#|          -the atomselections are processed to detect interacting residues .
#|          -selInfo array keys considered :-'<selId>,ligId', '<selId>,recId',
#|           _ '<selId>,ligSelTxt', and '<selId>,recSelTxt' .
#|            -'<selId>,frame' (optionally) ;
#|          -variables to be assigned (respectively): pdbIdL, pdbIdR,
#|           _ selTxtL, and selTxtR .
#|          -arguments that could override selInfo-stored values: 'pdbId*',
#|           _ 'frame', 'selTxt*' .
#|          -acceptable values :
#|            -an already registered selId for the global selInfo array .
#|            -"" :-the 'selId' arg will be ignored ;;
#|          -default value :-"" ;;
#|        -'pdbId', 'pdbFile', 'id', 'vmdId', 'molId' :
#|          -specs common for both lig and rec molecules .
#|          -can be interpreted as a 4-letter PDBID, a file name, or a VMD Id,
#|           _ depending on the 'source' argument .
#|          -this arg along with 'pdbIdL' and 'pdbIdR' may override the
#|           _ specs gotten from the selInfo array .
#|          -default value :-"top" ;;
#|        -'pdbIdL', 'pdbFileL', 'idL', 'vmdIdL', 'molIdL' :
#|          -same as arg 'pdbId' but only applies for the ligand molecule .
#|          -default value :-"" ;;
#|        -'pdbIdR', 'pdbFileR', 'idR', 'vmdIdR', 'molIdR' :
#|          -same as arg 'pdbId' but only applies for the receptor molecule .
#|          -default value :-"" ;;
#|        -'frame', 'frm', 'index' :
#|          -for .dlg files, it allows to specify the desired conformation .
#|          -defaukt value ;-0 ;;
#|        -'selTxtL' :
#|          -specifies the VMD's atom selection text for the ligand .
#|          -acceptable values :
#|            -valid selection text for the VMD's command 'atomselect' ;
#|          -default value :-"all" ;;
#|        -'selTxtR' :
#|          -specifies the VMD's atom selection text for the receptor atoms
#|          -overrides the atom selection text from the 'selInfo' array .
#|          -acceptable values :
#|            -valid selection text for the VMD's command 'atomselect' ;
#|          -default value :-"protein" ;;
#|        -'ligSel', 'atmSelL', 'atomSelLig' :
#|          -atom selection previously generated using the VMD's command
#|           _'atomselect' .
#|          -if specified atmSelL and atmSelR, will override other args
#|           _ specifying atom selection text such as selInfo or selTxt* .
#|          -default value :-"" ;;
#|        -'interactCAsel', 'iteractingCA', 'atmSelR', 'recSel',
#|         _ 'atomSelRec', 'atmSelCA' :
#|          -atom selection previously generated using the VMD's command
#|           _'atomselect' .
#|          -it must include only the alpha carbon atoms from the residues
#|           _ in the receptor in direct contact with the ligand .
#|          -if specified atmSelL and atmSelR, will override other args
#|           _ specifying atom selection text such as selInfo or selTxt* .
#|          -default value :-"" ;;
#|        -'workPath', 'workFolder', 'workDir',
#|         _ 'outPath', 'outFolder', 'outDir' :
#|          -directory path to be prepended to output files (before prefix) .
#|          -the path will be created in case it is not found ;
#|          -acceptable values :
#|            -a valid either relative or absolute directory path .
#|            -"", ".", "./" :
#|              -the current folder is used as workPath ;;
#|          -defautl value :-"" ;
#|        -'l_topFile', 'l_topFiles', 'topFile', 'topFiles',
#|         _ 'l_topologyFile', 'topologyFiles' :
#|          -list of topology files to be used by psfgen .
#|          -the file names may include either relative or absolute paths .
#|          -this list will be read before the top files found in the 'topDir' .
#|          -default value :-{} ;;
#|        -'topDir', 'topPath', 'topologyDir', 'topologyPath' :
#|          -path to the directory where topology files can be found .
#|          -acceptable values :
#|            -an absolute or relative dir path (ending with "/") .
#|            -a list of dir paths is also acceptable .
#|            -"" :-no directory is considered to search topology paths ;;
#|          -default value :-"" ;;
#|        -'l_topExt', 'l_topFileExt' :
#|          -list of known file extensions for CHARMM topology .
#|          -default value :-{top rtf str} ;;
#|        -'prefix', 'outPrefix', 'outputPrefix' :
#|          -prefix used for all pdb file generated .
#|          -the specified name may include abs or rel file paths .
#|          -acceptable values :
#|            -any string acceptable as part of a file name .
#|            -"auto" :- ;;
#|          -default value :-"auto" ;;
#|        -'pgnWrite', 'pgnScript', 'pgnFile', 'psfgenFile', 'psfgenScript' :
#|          -controls whether a psfgen is generated 'on-the-fly' or a psfgen
#|           _ config pgn file is written for later processing .
#|          -acceptable values :
#|            -0 :-the psfgen structures are created within this proc ;
#|            -1 :-a psfgen config pgn file is created instead ;;
#|          -default value :-1 ;;
#|        -'segPrefix', 'segmentPrefix', 'segLetter' :
#|          -one-to-three-letter name prefix for the created segments .
#|          -default value :-"F" ;;
#|        -'cutoff', 'distCutoff', 'distance', 'dist', 'r' :
#|          -maximum distance between lig and rec atoms to be considered as
#|           _ interaction .
#|          -default value :-4.5 ;;
#|        -'gapMax', 'maxGap', 'resGap', 'gap', 'gaps', 'gapSize', 'gapLen' :
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
#|          -default value :-2 ;;
#|        -'tailN', 'NTail' :
#|          -same as 'tails' arg but specific for the peptide N-terminus ;
#|        -'tailC', 'CTail' :
#|          -same as 'tails' arg but specific for the peptide C-terminus ;
#|        -'mutateGaps', 'gapsMutate', 'mutate' :
#|          -specifies whether gap residues should be mutated to Alanine .
#|          -***note: user-specified mutations may be considered in future
#|           _ versions .
#|          -acceptable values :-0 .-1 ;
#|          -default value :-1 ;
#|        -'keepGlycine', 'keepGapGly', 'glyGapKeep' :
#|          -specifies whether glycine rec-gap residues should be ignored or
#|           _ mutated according the 'mutateGaps' arg .
#|          -acceptable values :-0 .-1 ;
#|          -default value :-1 ;;
#|        -'keepProline', 'keepGapPro', 'proGapKeep' :
#|          -specifies whether proline rec-gap residues should be ignored or
#|           _ mutated according the 'mutateGaps' arg .
#|          -acceptable values :-0 .-1 ;
#|          -default value :-1 ;;
#|        -'keepCysteine, 'keepGapCys', 'cysGapKepp' :
#|          -specifies whether cysteine rec-gap residues should be ignored or
#|           _ mutated according the 'mutateGaps' arg .
#|          -acceptable values :-0 .-1 ;
#|          -default value :-0 ;;
#|        -'chargeTarget', 'targetCharge', 'mutateChargeTarget' :
#|          -***incomplete implementation*** .
#|          -acceptable values :
#|            -"0", "zero", "neutral" :
#|              -a neutral charge for the whole system is sought .
#|              -two options will be implemented in future versions :
#|                -all residues neutral .
#|                -other strategies ;;
#|            -"" :-the chargeTarget arg will be ignored ;;
#|          -default value :-"" ;;
#|        -'bll', 'baseLogLevel', 'baseLogLvl' :
#|          -base log level used for log messages in the proc .
#|          -see library logLib .
#|          -defaul value :-1 ;;
#|        -'resAtm', 'pivotTxt', 'pivotSelTxt', 'uniqueResAtom' :
#|          -selection text specifying one atom appearing once on each residue .
#|          -default value :-"name CA" ;;
#|        -'nTerPatch', 'firstPatch' :
#|          -N-terminal patch (PRES) name (psfgen keyword 'first') .
#|          -default value :-"NNEU" ;;
#|        -'cTerPatch' 'lastPatch' :
#|          -C-terminal patch (PRES) name (psfgen keyword 'last') .
#|          -default value :-"CNEU" ;;
#|        -'nTerGlyPatch', 'nTerPatchGly', 'firstPatchGly', 'firstGlyPatch' :
#|          -N-terminal patch (PRES) name for Gly (psfgen keyword 'first') .
#|          -default value :-"NGNE" ;;
#|        -'nTerProPatch', 'nTerPatchPro' 'firstPatchPro', 'firstProPatch' :
#|          -N-terminal patch (PRES) name for Gly (psfgen keyword 'first') .
#|          -default value :-"ACP" ;;
#|        - ;;
#|    -notes :
#|      -'pdbAliasCad' fixed .
#|      -topology file/dir specification mandatory (no default values) ;;

proc lr_trimComplex {complexType args} {
# global variables
  global selInfo
# configuring logLib to manage log messages
  namespace import ::logLib::*
  set procName [lindex [info level 0] 0]
  set_logName $procName
  set_logVersion "002"
  set bll 1   ;# base log level
# interpreting logLib arguments
  set args_rest [eval arg_interpreter $args]
  logMsg "+++ Trim of ligand-receptor complexes +++" $bll
  logMsg " command line: $procName $complexType $args" $bll
# default value for variables
  set src "id"
  set pdbId "top"
  set pdbIdL ""
  set pdbIdR ""
  set frame 0
  set selId ""
  set selTxtL "all"
  set selTxtR "protein"
  set ligSel ""
  set interactCASel ""
  set workPath ""
  set l_topFile {}
  set topDir ""
  set l_topExt {top rtf str}
  set prefix "auto"
  set pgnWrite 1
  set segPrefix "F"
  set cutoff 4.5
  set gapMax 3
  set tails 2
  set tailN $tails
  set tailC $tails
  set mutateGaps 1
  set keepGlycine 1
  set keepProline 1
  set keepCysteine 0
  set chargeTarget ""
  set resAtm "name CA"
  set nTerPatch "NNEU"
  set cTerPatch "CNEU"
  set nTerGlyPatch "NGNE"
  set nTerProPatch "ACP"   ;# acetylated N-ter proline (neutral)
  set pdbAliasCad "pdbalias residue HIS HSD\npdbalias atom ILE CD1 CD\npdbalias residue HOH TP3M\npdbalias atom TP3M O OH2\n"
  set args_last {}
# read input command-line arguments (excluding logLib arguments)
  if {[expr {[llength $args_rest]%2}] == 0} {
    foreach {arg val} $args_rest {
      switch [string tolower $arg] {
        "src" - "source" - "input" - "inputsource" {set src $val}
        "pdbid" - "pdbfile" - "id" - "vmdid" - "molid" {set pdbId $val}
        "pdbidl" - "pdbfilel" - "idl" - "vmdidl" - "molidl" {set pdbIdL $val}
        "pdbidr" - "pdbfiler" - "idr" - "vmdidr" - "molidr" {set pdbIdR $val}
        "frame" - "frm" - "index" {set frame $val}
        "selid" - "selinfo" - "selinfoid" - "selinfokey" {set selId $val}
        "seltxtl" {set selTxtL $val}
        "seltxtr" {set selTxtR $val}
        "ligsel" - "atmsell" - "atomsellig" {set ligSel $val}
        "interactcasel" - "atmselr" - "interactingca" - "recsel" \
          - "atomselrec" - "atmselca" {set interactCASel $val}
        "workpath" - "workfolder" - "workdir" - "outpath" - "outfolder" \
          - "outdir" {set workPath $val}
        "l_topfile" - "l_topfiles" - "topfile" - "topfiles" \
          - "l_topologyfile" - "topologyfiles" {set l_topFile $val}
        "topdir" - "toppar" - "topologydir" - "topologypath" {set topDir $val}
        "l_topExt" - "l_topFileExt" {set l_topExt $val}
        "prefix" - "outprefix" - "outputprefix" {set prefix $val}
        "pgnwrite" - "pgnscript" - "pgnfile" - "psfgenfile" - "psfgenscript" \
          {set pgnWrite $val}
        "segprefix" - "segmentprefix" - "segletter" {set segPrefix $val}
        "cutoff" - "distcutoff" - "distance" - "dist" - "r" {set cutoff $val}
        "gapmax" - "maxgap" - "resgap" - "gap" - "gaps" - "gapsize" - \
          "gaplen" - "lengap" - "gaplength" - "lengthgap" {set gapMax $val}
        "tails" - "tailslength" - "lengthtails" - "tailslen" - "lentails" \
          {set tails $val}
        "tailn" - "ntail" {set tailN $val}
        "tailc" - "ctail" {set tailC $val}
        "mutategaps" - "gapsmutate" - "mutate" {set mutateGaps $val}
        "keepglycine" - "keepgapgly" - "glygapkeep" {set keepGlycine $val}
        "keepproline" - "keepgappro" - "progapkeep" {set keepProline $val}
        "keepcysteine" - "keepgapcys" - "cysgapkeep" {set keepCysteine $val}
        "chargetarget" - "targetcharge" - "mutatechargetarget" {}
        "bll" - "baseloglevel" - "baseloglvl" {set bll $val}
        "resatm" - "pivottxt" -  "pivotseltxt" - "uniqueresatom" \
          {set resAtm $val}
        "nterpatch" - "firstpatch" {set nTerPatch $val}
        "cterpatch" - "lastpatch" {set cTerPatch $val}
        "nterglypatch" - "firstpatchgly" {set nTerGlyPatch $val}
        "nterpropatch" - "firstpatchpro" {set nTerProPatch $val}
        default {
          logMsg "$procName: argument (value) unknown: $arg ($val)" $bll
          set args_last [concat ${args_last} $arg $val]
          }
        }
      }
    logMsg "$procName: unused arguments: $args_last"
  } else {
    logMsg "$procName: Odd number of arguments: $args_rest" $bll
    return ""
    }
# updating log levels
  set ll1 $bll
  set ll2 [expr {$bll+1}]
  set ll3 [expr {$bll+2}]
  logMsg " using base log level: $bll" $ll2
# process input variable values
# check complexType, src, pdbId* options and load molecules
  switch [string tolower $complexType] {
    "pdb" {   ;# single PDB file with both the lig and the rec
      set singleMol 1
      switch [string tolower $src] {
        "download" - "pdbid" - "pdbload" - "pdb" {
          set id [mol pdbload $pdbId]
          set idL $id; set idR $id
          logMsg " PDB model downloaded into id $id: $pdbId" $ll2
          }
        "file" - "loadfile" - "localfile" - "localpdb" {
          set id [mol new $pdbId type pdb waitfor all]
          set idL $id; set idR $id
          logMsg " loaded local file into id $id: $pdbId" $ll2
          animate goto start
          set name [string range $pdbId \
            [expr {[string last / $pdbId] + 1}] \
            [expr {[string last ".pdb" $pdbId] - 1}]]
          mol rename $id $name
          logMsg " mol renamed to: $name" $ll3
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
          if {$selTxtL == "all"} {
            set selTxtL $selInfo($selId,ligSelTxt)
            }
          if {$selTxtR == "protein"} {
            set selTxtR $selInfo($selId,recSelTxt)
            }
          logMsg " setting selTxtL: $selTxtL" $ll3
          logMsg " setting selTxtR: $selTxtR" $ll3
          if {[info exists selInfo($selId,frame)]} {
            set frame $selInfo($selId,frame)
            logMsg " setting frame: $frame" $ll3
            }
          if {$selInfo($selId,ligId) == $selInfo($selId,recId)} {
            set singleMol 1
            set id $selInfo($selId,ligId)
            set idL $id; set idR $id
            logMsg " setting id (same for both lig and rec): $id" $ll3
          } else {
            set singleMol 0
            set idL $selInfo($selId,ligId)
            set idR $selInfo($selId,recId)
            logMsg " setting idL: $idL" $ll3
            logMsg " setting idR: $idR" $ll3
            }
        } else {
          logMsg "$procName: selInfo missing data for selId: $selId"
          logMsg ", required array keys ($selId,<key>): ligId, recId, ligSelTxt, recSelTxt, and title" $ll1
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
      lr_trimComplex_help [lindex $args 0]
      return ""
      }
    default {
      logMsg "Unknown option for 'complexType' arg: $src" $ll1
      return ""
      }
    }
# create atom selections
  if {$singleMol} {
    if {($ligSel == "") || ($interactCASel == "")} {
      logMsg " using interaction cutoff value: $cutoff" $ll2
      logMsg " using selTxtL: $selTxtL" $ll2
      logMsg " using selTxtR: $selTxtR" $ll2
      logMsg " using resAtm: $resAtm" $ll2
      set ligSel [atomselect $id "$selTxtL" frame $frame]
      set interactCASel [atomselect $id "(same residue as protein within $cutoff of $selTxtL) and $selTxtR and $resAtm"]
    } else {
      logMsg " atom selections previously defined:" $ill2
      logMsg "  lig: [$ligSel text] " $ill2
      logMsg "  rec: [$interactCASel text]" $ll2  
      }
    logMsg "rec residues selected: \
            [$interactCASel get {chain resname resid}]" $ll2
  } else {
    # not implemented yet
    logMsg "lig and rec from different molecules not yet implemeted" $ll1
    return ""
    }
# check workPath, prefix, and pgnWrite arguments
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
  if {$prefix == "auto"} {
    set prefix "[molinfo $idR get name]-[lindex [$ligSel get resname] 0]_${cutoff}"
    logMsg " generated prefix: $prefix" $ll3
    }
  logMsg " setting prefix for output files: $prefix" $ll2
  if {$pgnWrite} {
    set pgnOut [open "${workPath}${prefix}.pgn.tcl" w]
    logMsg "" $ll2
  } else {
    logMsg " pgn script for psfgen will not be written" $ll2
    }
# read topology files
  if {($l_topFile == {}) && ($topDir == "")} {
    logMsg "'l_topFile' and/or 'topDir' arguments must be specified" $ll1
    return ""
    }
  set topCad ""
  if {[llength $l_topFile] >= 1} {
    logMsg "reading topology files from l_topFile: $l_topFile" $ll2
    foreach topFile $l_topFile {
      if {[file exists $topFile]} {
        logMsg "topology file found: $topFile" $ll3
        set topCad "${topCad}topology $topFile\n"
      } else {
        logMsg "topology file not found: $topFile" $ll2
        }
      }
    }
  if {$topDir != ""} {
    foreach dir $topDir {
      logMsg "searching for topology files (ext: $l_topExt) in dir: $dir" $ll2
      foreach ext $l_topExt {
        if {[catch {lsort [eval glob "${dir}*.${ext}"]} topDirFiles]} {
          set topDirFiles ""
          logMsg "no .${ext} files found in dir ${dir}" $ll3
          }
        foreach topFile $topDirFiles {
          logMsg "topology file found: $topFile" $ll3
          set topCad "${topCad}topology $topFile\n"
          }
        }
      }
    }
  if {$topCad == ""} {
    logMsg "no topology files were found" $ll1
    return ""
    }
# detect chain fragments
  logMsg " using parameters for fragment building:" $ll2
  logMsg "  gapMax: $gapMax  tailN: $tailN  tailC: $tailC" $ll2
  logMsg "  mutateGaps: $mutateGaps  " $ll2
  logMsg "  keepProline: $keepProline  keepCysteine: $keepCysteine" $ll2
  logMsg "  chargeTarget: $chargeTarget" $ll3
  logMsg "  NTerPatch: $nTerPatch  CTerPatch: $cTerPatch" $ll2
  logMsg "  NTerGlyPatch: $nTerGlyPatch  NTerProPatch: $nTerProPatch" $ll3
# classify resids by chain
  if {$singleMol} {
    foreach chain [$interactCASel get chain] rid [$interactCASel get resid] {
      if {[info exists resids($chain)]} {
        lappend resids($chain) $rid
      } else {
        set resids($chain) $rid
        }
      }
    logMsg "chains in contact with the ligand(s): [array names resids]" $ll2
# header for pgn file
    if {$pgnWrite} {
      puts $pgnOut "# psfgen config pgn file created by: [get_logName_version]"
      puts $pgnOut ""
      puts $pgnOut $topCad
      puts $pgnOut $pdbAliasCad
    } else {
      package require psfgen
      eval $topCad
      eval $pdbAliasCad
      }
    logMsg " psfgen header: topCad:" $ll2
    logMsg $topCad $ll2
    logMsg " psfgen header: pdbAliasCad:" $ll2
    logMsg $pdbAliasCad $ll2
    # process each chain
    set ll_frag_rid {}
    foreach chain [array names resids] {
      set prevrid 0   ;# last resid assigned to a fragment (within a chain)
      set l_rid {}   ;# list of resids for the current fragment
      set l_gapTail {}   ;# list of resids added to the current gap
      set nFrag 1
      set r 0
      while {r < [llength $resids($chain)]} {
        set rid [lindex $resids($chain) $r]
        set seg "${chain}${nFrag}"
        # filling a gap or adding an N-tail
        if {[llength $l_rid] == 0} {   ;# first interacting resid
          # adds an N-tail if there are residues available
          set i [expr {$rid-$tailN}]
          while {$i < $rid} {
            if {$nFrag == 1} {   ;# starting segment without N-tail
              set i $rid   ;# makes this loop run only once
              incr r
              }
            set gapId [atomselect $idR "chain $chain and resid $i and $resAtm"]
            if {[$gapId num] == 1} {
              # start psfgen segment and add first residue
              logMsg "starting segment ${seg}" $ll3
              logMsg "added resid $i as N-tail" $ll3
              if {$pgnWrite} {
                puts $pgnOut "segment ${seg} \{"
                puts $pgnOut "  residue $i [$gapId get resname] $chain"
              } else {
                eval "segment ${seg} \{"
                eval "  residue $i [$gapId get resname] $chain"
                }
              lappend l_rid $i
              # mutate fist residue?
              switch [$gapId get name] {
                "GLY" {
                  if {$keepGlycine} {
                    logMsg " first patch for segment ${seg}: $nTerGlyPatch" $ll3
                    if {$pgnWrite} {puts $pgnOut "  first $nTerGlyPatch"} else {eval "first $nTerGlyPatch"}
                    }
                  }
                "PRO" {
                  if {$keepProline} {
                    logMsg " first patch for segment ${seg}: $nTerProPatch" $ll3
                    if {$pgnWrite} {puts $pgnOut "  first $nTerProPatch"} else {eval "first $nTerProPatch"}
                    }
                  }
                "CYS" {
                  if {$keepCysteine} {
                    logMsg " first patch for segment ${seg}: $nTerPatch" $ll3
                    if {$pgnWrite} {puts $pgnOut "  first $nTerPatch"} else {eval "first $nTerPatch"}
                    }
                  }
                default {
                  logMsg " first patch for segment ${seg}: $nTerPatch" $ll3
                  if {$mutateGaps} {logMsg " residue $i mutated to ALA" $ll3}
                  if {$pgnWrite} {
                    puts $pgnOut "  first $nTerPatch"
                    if {$mutateGaps} {puts $pgnOut "  mutate $i ALA"}
                  } else {
                    eval "  first $nTerPatch"
                    if {$mutateGaps} {eval "  mutate $i ALA"}
                    }
                  }
                }
              $gapId delete
              break
            } else {
              logMsg "unable to add resid $i as N-tail" $ll3
              }
            $gapId delete
            incr i
            }   ;# while
        } else {   ;# next interacting resid within a chain/fragment
# *** check adding first interaction sphere resids***
          set prevrid [lindex l_rid end]
          if {$rid == [expr {$prevrid+1}]} {  ;# consecutive resid
            logMsg "adding residue to segment ${seg}: [[atomselect $idR "chain $chain and resid $rid and $resAtm"] get resname] $rid" $ll3
            if {$pgnWrite} {
              puts $pgnOut "  residue $rid [$gapId get resname] $chain"
            } else {
              eval "residue $rid [$gapId get resname] $chain"
              }
            lappend l_rid $rid
            incr r
          } else {   ;# non-consecutive resid
            if {$prevrid >= [expr {$rid-$gapMax-1}]} {   ;# add gap resid 
              for {set i [expr {$prevrid+1}]} {$i < $rid} {incr i} {
                set gapId [atomselect $idR "chain $chain and resid $i and $resAtm"]
                if {[$gapId num] == 1} {
                  logMsg "adding residue [$gapId get resname] $i to segment ${chain}${nFrag}" $ll3
                  if {$pgnWrite} {
                    puts $pgnOut "  residue $i [$gapId get resname] $chain"
                  } else {
                    eval "residue $i [$gapId get resname] $chain"
                  }
                  lappend l_rid $i
                  # mutate gap resid?
                  switch [$gapId get resname] {
                    "GLY" {
                      if {!$keepGlycine} {
                        logMsg "GLY $i mutated to ALA in seg ${chain}${nFrag}" $ll3
                        if {$pgnWrite} {puts $pgnOut "  mutate $i ALA"} else {eval "  mutate $i ALA"}}}
                    "PRO" {
                      if {!$keepProline} {
                        logMsg "PRO $i mutated to ALA in seg ${chain}${nFrag}" $ll3
                        if {$pgnWrite} {puts $pgnOut "  mutate $i ALA"} else {eval "  mutate $i ALA"}}}
                    "CYS" {
                      if {!$keepCysteine} {
                        logMsg "CYS $i mutated to ALA in seg ${chain}${nFrag}" $ll3
                        if {$pgnWrite} {puts $pgnOut "  mutate $i ALA"} else {eval "  mutate $i ALA"}}}
                    default {
                      logMsg "[$gapId get resname] $i mutated to ALA in seg ${chain}${nFrag}" $ll3
                      
                      if {$pgnWrite} {puts $pgnOut "  mutate $i ALA"} else {eval "  mutate $i ALA"}
                      }
                    }
                } else {
                  logMsg "missing residue in chain $chain: $i" $ll1
                  return ""
                  }
                }
            } else {   ;# new fragment should be added
              logMsg " added last (cTer) patch: $cTerPatch" $ll2
              logMsg " closing segment: ${chain}${nFrag}" $ll2
              if {$pgnWrite} {
                puts $pgnOut "  last $cTerPatch"
                puts $pgnOut "  \}"
              } else {
                eval "last $cTerPatch"
                eval "\}"
                }
              logMsg " resids for complete segment ${seg}: $l_rid" $ll3
              lappend ll_frag_rid $ll_rid
              logMsg " current ll_frag_rid: ${ll_frag_rid}" $ll3
              set l_rid {}
              incr nFrag
              }
            }
          }
        }   ;# while
      # adding tailC
      
      }
  } else {
    logMsg "chain fragment detection not yet implemented for two molecs " $ll1
    return ""
    }
  close $pgnOut

return ""

# report working values for variables

set res $atmSelR
# process fragments
  set nFrag 1
  set prevrid 0
  set prefName $prefix
  set nResFinal 1
  set l_gaps {}
  set l_rid {}
  set l_fragGaps {}

foreach rid [$res get resid] {
# search for consecutive residues to form a single fragment
# note: will not work if $res starts in 1 or 2
  if {$rid == [expr {$prevrid + 1}]} {
    lappend l_rid $rid
    incr nResFinal
  } elseif {$rid == [expr {$prevrid + 2}]} {
    lappend l_rid [expr {$rid - 1}]
    lappend l_rid $rid
    puts "residue gap added: [expr {$rid - 1}]"
    lappend l_gaps [expr {$rid - 1}]
    lappend l_fragGaps [expr {$rid - 1}]
    incr nResFinal 2
  } elseif {$rid == [expr {$prevrid + 3}]} {
    lappend l_rid [expr {$rid - 2}]
    lappend l_rid [expr {$rid - 1}]
    lappend l_rid $rid
    puts "residue gap added: [expr {$rid - 2}]"
    puts "residue gap added: [expr {$rid - 1}]"
    lappend l_gaps [expr {$rid - 2}]
    lappend l_gaps [expr {$rid - 1}]
    lappend l_fragGaps [expr {$rid - 2}]
    lappend l_fragGaps [expr {$rid - 1}]
    incr nResFinal 3
  } elseif {$rid == [expr {$prevrid + 4}]} {
    lappend l_rid [expr {$rid - 3}]
    lappend l_rid [expr {$rid - 2}]
    lappend l_rid [expr {$rid - 1}]
    lappend l_rid $rid
    puts "residue gap added: [expr {$rid - 3}]"
    puts "residue gap added: [expr {$rid - 2}]"
    puts "residue gap added: [expr {$rid - 1}]"
    lappend l_gaps [expr {$rid - 3}]
    lappend l_gaps [expr {$rid - 2}]
    lappend l_gaps [expr {$rid - 1}]
    lappend l_fragGaps [expr {$rid - 3}]
    lappend l_fragGaps [expr {$rid - 2}]
    lappend l_fragGaps [expr {$rid - 1}]
    incr nResFinal 4
    incr nResFinal 4
  } else {
    if {[llength $l_rid] == 0} {
      set l_rid $rid
      incr nResFinal
    } else {
      set rsel [atomselect top "resid $l_rid and chain A and not altloc B"]
      $rsel writepdb ${prefName}_frag${nFrag}.pdb
      puts "writing fragment $nFrag, residues $l_rid"
      $rsel delete
      lappend ll_gaps $l_fragGaps
      incr nFrag
      set l_rid $rid
      incr nResFinal
  set l_fragGaps {}
      }
    }
  set prevrid $rid
  }
  if {[llength $l_rid] != 0} {
    set rsel [atomselect top "resid $l_rid and chain A and not altloc B"]
    $rsel writepdb ${prefName}_frag${nFrag}.pdb
    puts "writing fragment $nFrag, residues $l_rid"
    }
  puts "gap residues added: $l_gaps"
  puts "Total number of residues: $nResFinal"
puts "ll_gaps $ll_gaps"

package require psfgen

topology ../complexes/toppar/top_all36_prot.rtf
topology ../complexes/toppar/top_all36_cgenff.rtf
topology ../complexes/toppar/OLCprot.str
topology ../complexes/toppar/91Qprot.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
pdbalias residue HOH TP3M
pdbalias atom TP3M O OH2



  set nfrag 0
  foreach l_frag $ll_gaps {
    incr nfrag
    set ngap 0
    
    foreach gap $l_frag {
      incr ngap
      set gapSel [atomselect top "$selTxtR and resid $gap"]
      puts "gapSelTxt [$gapSel text] "
      if {$ngap == 1} {
        if {[$gapSel get resname] == "GLY"} {
        } elseif {[$gapSel get resname] == "PRO"} {
        } else {
          }
        }
      }
    }

  
  }   ;# proc lr_trimComplex

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



#|- ;





