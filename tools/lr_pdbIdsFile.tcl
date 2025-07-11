
#|-lr_pdbIdsFile {l_pdbId {src "download"} args} :
#|  -detect ligand-receptor complexes in lists of PDB files .
#|  -this procedure is part of the https://github.com/carloszep/lrcd-vmd.tcl
#|   _ library .
#|  -creates a pdbIdsFile used by proc lrcdVec_writeDB to store LRC vectors .
#|  -pdbIdsFile format (4 columns) :
#|     -(1) 4-letter PDBID .
#|     -(2) 1-letter chain identifier (i.e. A) .
#|     -(3) 3-letter PDB residue name for the ligand .
#|     -(4) 1- to 4-letter residue seq. number (resid) ;
#|  -for each PDBID and chain, all non-protein residues are browsed
#|   _ to detect "native" LR complexes to populate a list of ligands .
#|  -optionally a folder tree is created for cleaned structures of proteins,
#|   _ ligands, and complexes (for docking VS calculations) .
#|  -the argument 'workPath' must be specified to create the folder tree .
#|  -folder tree and files created (for automated autodock4 calculations) :
#|    -'<workPath>' :
#|      -'lig/' :
#|        -'pdb/' :
#|          -'refs/' :
#|            -'<resName><resId>-<pdbId><chain>.pdb' ;;
#|        -'nopt/' :
#|          -'noh/' :
#|            -'pdb/' :
#|              -'<resName><resId>-<pdbId><chain>.pdb' ;;;;
#|      -'rec/' :
#|        -'pdb/' :
#|          -'full/' :
#|            -'<pdbId>.pdb' ;
#|          -'rcsb.org/' :
#|            -'<pdbId>.pdb' ;
#|          -'recs/' :
#|            -'<pdbId>.pdb' ;
#|          -'<pdbId>.pdb' ;;
#|      -'grid/' :
#|        -'<pdbId><chain>-<resName><resId>/' ;
#|      -'dock/' :
#|        -'<pdbId><chain>-<resName><resId>/' ;;;
#|  -arguments :
#|    -l_pdbId :
#|      -list of PDBIDs to be processed .
#|      -the source of the PDB files depends on the value of the 'src' arg .
#|    -src :
#|      -specify the source of the PDB files to be processed .
#|      -acceptable values :
#|        -"download", "pdbLoad", "internet" :
#|          -the l_pdbid arg must be a list of PDBIDs to be downloaded
#|           _ from the RCSB PDB web site (uses 'mol pdbload') ;
#|        -"file", "loadFile", "pdbFile" :
#|          -the l_pdbid arg must be a list of PDB file names to be loaded .
#|          -PDB format is assumed .
#|          -it is recommended to avoid including paths to files, instead
#|           _ use the argument 'pdbPath' ;
#|        -"molId", "idl", "id", "ids", "vmdId", "loaded" :
#|          -the l_pdbid arg must be a list of preloaded VMD mol Ids ;;
#|      -default value :-'download' ;;
#|    -args (variable arguments, case insensitive) :
#|      -'pdbIdsFile', 'out', 'output', '-o' :
#|        -name of the file for writing the list of Ids of ligands .
#|        -if the file already exists the output will be appended .
#|        -default value :-"pdbId-ligDB.txt" ;;
#|      -'workPath', 'workFolder', 'outTree', 'dirTree' :
#|        -path prepended to all output files and folders .
#|        -the path must end with the '/' character .
#|        -acceptable values :
#|          -"none" :
#|            -no output folders for cleaned structures are created ;
#|          -"", ".", "./" :
#|            -the current folder is used as workPath ;
#|          -an absolute or relative (linux) folder path finishing with '/' ;
#|        -default value :
#|          -"none" ;;
#|      -'chain', 'chains', 'l_chain' :
#|        -list of chain identifiers to be considered .
#|        -acceptable values :
#|          -list of chain identifiers (e.g. '{A B P}') :
#|            -overrides automatically browsing all chains ;
#|          -'all', {}, "" :
#|            -all chains contained in the PDB are considered ;
#|        -default value :-'all' ;;
#|      -'resName', 'resNames',  'l_resName' :
#|        -list of non-protein residue names to be considered as ligands .
#|        -resnames for water molecules are ignored (e.g. WAT, HOH) .
#|        -acceptable values :
#|          -list of 3-letter residue names for ligands (e.g. '{ACD NAG HEM}') :
#|            -overrides automatically detecting ligand residue names ;
#|          -'all', {}, "" :
#|            -all non-protein resnames (but waters) are considered ;;
#|        -default value :-'all' ;
#|      -'resId', 'resIds', 'l_resId' :
#|        -list of residue seq. numbers for ligands to be considered .
#|        -acceptable values :
#|          -list of 1- to 4-letter resid numbers for ligands (e.g. {600 700}) .
#|          -'all' :
#|            -all non-protein resids are considered ;;;
#|      -'ll_ligRec', 'ligRecs' :
#|        -specification of non-protein ligands to be considered as part of
#|         _ the receptor structure .
#|        -formatted as a list of lists .
#|        -each sublist contains the specification of a single ligand .
#|        -the specified ligands will be excluded as ligands in complexes .
#|        -list format :
#|          -{{<pdbId> <chain> <resName> <resId>} ...} :
#|            -<pdbId> may take the value '".*"' to specify all the pdbIds
#|             _ but it is not a real regular expression .
#|            -<chain> and <resName> may contain regular expressions :
#|              -examples: '\".*\"', '\"[ABCD]\"', ',\"^C.*\"' ;
#|            -<resId> may be either a specific value or a range
#|             _ (e.g. '"1 to 100"', '192') .
#|            -examples :
#|              -{{".*" \"[ABCD]\" \"C[AB][0-9]+\" "1 to 100"}} .
#|              -{{6cox A HEM 682} {3pgh "A B C D" FLP 701}} ;;;;
#|      -'ll_ligExclude', 'ligExcludes' :
#|        -list of lists with the specification of ligands to be excluded .
#|        -each sublist contains the specification of a single ligand .
#|        -list format :-see 'll_ligRec' argument above ;;
#|      -'chainRecExclude', 'recChainExclude', 'recExclude' :
#|        -chain dentifiers to be excluded from the receptor molecule .
#|        -affects only the 'rec/pdb/<pdbId>.pdb' file .
#|        -acceptable values :
#|          -"none", {}, "" :-no chains are excluded ;
#|          -list of chain identifiers (e.g. "A", "B C E P") ;
#|        -default value :-{} ;;
#|      -'pdbPath', 'molPath' :
#|        -common folder path where the PDB files are loaded from .
#|        -used only with the argument 'source' set to "file" :
#|          -the l_pdbId argument would be expected to contain only file names
#|           _ without a path ;
#|        -the file path must finish with a '/' character .
#|        -example: "/path/to/files/"
#|        -default value :-"" ;;
#|      -'minLigSize', 'minLigNumAtm' :
#|        -the ligands detected will be excluded if has less atoms than the
#|         _ value of this parameter .
#|        -acceptable values :
#|          -positive integer ;
#|        -default value :
#|          -7 ;;
#|      -'loSt', 'channelId', 'log' :
#|        -output stream (channel) for log messages .
#|        -acceptable values :
#|          -a channel Id refering to an already opened output text file .
#|          -"none":
#|            -will avoid any output ;;
#|        -default value :
#|          -stdout ;;;;
#|  -notes :
#|    -pending to fix adding rec ligads through ll_ligRec i.e. COH in 5kir .
#|    -if no ligand is found in a chain, a grid is created centered
#|     _ in the chain's COM .
#|    -receptor pdb files are duplicated in the rec/pdb and rec/pdb/recs dirs .
#|    -from version 1.0.6 the resid of a ligand is included for a total of
#|     _ four columns in the pdbIdFile .
#|    -water molecules are excluded by default .
#|    -from v111, the parameter minLigNumAtm and the list exclLigName will
#|     _ define the filtering of ligands .
#|    -in case a ligand is interacting with more than one chain, the receptor
#|     _ moiety should consider a list of chain ids :
#|      -to be considered when writing receptor structures for docking ;
#|    -the vmd 'name' field for the pdbs must be the same as the PDBID .
#|    -exclLigName stores a list residue names to be ignored :
#|      -current list: "HOH" "WAT" "PEG" "NAG" "FUC" "FUL" ;
#|    -for the moment, no frame specification in each pdb is handled :
#|      -after loading the molecule, the first frame is seeked ;;;
proc lr_pdbIdsFile {l_pdbId {src "download"} args} {
# global variables ...
# default values for variables and arguments
  set procName [lindex [info level 0] 0]
  set pdbIdsFile "pdbId-ligDB.txt"
  set l_chainUsr {}
  set l_resNameUsr {}
  set l_resIdUsr {}
  set ll_ligRec {}
  set ll_ligExclude {}
  set recExcl {}
  set pdbPath ""
  set workPath "none"
  set loSt stdout
  set out 1
  set minLigNumAtm 7   ;#minimum num of atoms to consider a residue as a ligand
  set exclLigName [list "HOH" "WAT" "PEG" "NAG" "FUC" "FUL" "MAN" "BMA" "PGE" "GLC"]
  set args_rest {}   ;# list of unrecognized arguments; may be reused
# decode variable arguments
  if {[expr {[llength $args]%2}] == 0} {   ;# even or 0 optional arguments
    foreach {arg val} $args {
      switch [string tolower $arg] {
        "pdbidsfile" - "output" - "out" - "-o" {set pdbIdsFile $val}
        "workpath" - "workfolder" - "outtree" - "dirtree" {set workPath $val}
        "chain" - "chains" - "l_chain" {set l_chainUsr $val}
        "resname" - "resnames" - "l_resname" {set l_resNameUsr $val}
        "resid" - "resids" - "l_resid" {set l_resIdUsr $val}
        "ll_ligrec" - "ligrecs" {set ll_ligRec $val}
        "ll_ligexclude" - "ligexcludes" {set ll_ligExclude $val}
        "recchainexclude" - "chainrecexclude" - "recexclude" {set recExcl $val}
        "pdbpath" - "molpath" {set pdbPath $val}
        "minligsize" - "minlignumatm" {set minLigNumAtm $val}
        "loSt" - "lost" - "channelId" - "channelid" - "log" {
          set loSt $val
          if {$loSt == "none"} {set out 0}
          }
        default {
          if $out {puts $loSt "$procName: argument unkown: $arg ($val)"}
          set args_rest [concat ${args_rest} $arg $val]
          }
        }
      }
    } else {   ;# odd number of arguments
      if $out {puts $loSt "$procName: Odd number of arguments: $args"}
      return ""
      }
# report user-provided parameters (arguments) to the log output ...
  if $out {
    puts $loSt "\n+++ Autodetect ligand-receptor complexes in PDBs +++"
    puts $loSt " command line: $procName {$l_pdbId} $src $args"
    puts $loSt " list of pdbIds (source: $src): ${l_pdbId}"
    puts $loSt " output pdbIdsFile: $pdbIdsFile"
    puts $loSt " user list of non-prot res in the receptor: $ll_ligRec"
    puts $loSt " user list of non-prot res to be excluded: $ll_ligExclude"
    puts $loSt " chains excluded from the output receptor file: $recExcl"
    puts $loSt " minimum number of atoms in ligands: $minLigNumAtm"
    puts $loSt " default list of resNames excluded as ligands: $exclLigName"
    puts $loSt " pdbPath for input files: $pdbPath"
    puts $loSt " output workPath: $workPath"
    if {$args_rest != {}} {puts $loSt " unused arguments: $args_rest"}
    }
# loading pdb molecules, and creating list of ids
  set l_id {}
  switch [string tolower $src] {
    "download" - "pdbload" - "internet" {
      foreach pdbId ${l_pdbId} {
        lappend l_id [mol pdbload $pdbId]
        animate goto start
        }
      }
    "file" - "loadfile" {
      foreach pdbName ${l_pdbId} {
        set pdbFile ${pdbPath}${pdbName}
        set id [mol new $pdbFile type pdb waitfor all]
        animate goto start
        set name [string range $pdbFile \
          [expr {[string last / $pdbFile] + 1}] \
          [expr {[string last ".pdb" $pdbFile] - 1}]]
        mol rename $id $name
        lappend l_id $id
        }
      }
    "molid" - "idl" - "id" - "ids" - "vmdId" - "loaded" {set l_id ${l_pdbId}}
    default {
      if $out {puts $loSt "$procName: unknown option in 'source' arg."}
      return ""
      }
    }
# creating output files and folders
  set outIds [open $pdbIdsFile a]
  if {$workPath == "none"} {set saveDir 0} else {set saveDir 1}
  if {($workPath == ".") || ($workPath == "")} {set workPath "./"}
  if {[string index $workPath end] != "/"} {set workPath "$workPath/"}
  if $saveDir {
    if {![file exist $workPath]} {exec mkdir $workPath}
  } else {
    if $out {
      puts $loSt "  Warning: no 'workPath' specified, no output directory tree for VS will be created."
      }
    }
# processing each mol id
  foreach id ${l_id} {
    set l_chain $l_chainUsr
    set pdbId [molinfo $id get name]
    if $saveDir {
      exec mkdir -p "${workPath}"
      exec mkdir -p "${workPath}lig/pdb/refs/"
      exec mkdir -p "${workPath}lig/pdbqt/"
      exec mkdir -p "${workPath}lig/nopt/noh/pdb/"
      exec echo "Use genligpdbqt.sh to populate from ../pdb/*.pdb files." > "${workPath}lig/pdbqt/README.txt"
      exec mkdir -p "${workPath}rec/pdb/recs/"
      exec mkdir -p "${workPath}rec/pdbqt/"
      exec mkdir -p "${workPath}rec/pdb/full/"
      [atomselect $id "all"] writepdb "${workPath}rec/pdb/full/${pdbId}.pdb"
      exec echo "Use genrecpdbqt.sh to populate from ../pdb/*.pdb files." > "${workPath}rec/pdbqt/README.txt"
#      catch {exec wget http://files.rcsb.org/download/${pdbId}.pdb} fid
      catch {exec curl -o ${pdbId} "http://files.rcsb.org/view/${pdbId}.pdb"} fid
      if {[file exist ${pdbId}.pdb]} {
        exec mkdir -p "${workPath}rec/pdb/rcsb.org/"
        exec mv "${pdbId}.pdb" "${workPath}rec/pdb/rcsb.org/"
      } else {
        puts $loSt "  Could not download file ${pdbId}.pdb from rcsb.org."
        }
      }
    set selTxtLigRec "(protein or nucleic)"
    if {(${l_chain} == "all") || (${l_chain} == {})} {
# look for available chains
      set l_chain {}
      foreach indAt [[atomselect $id "all"] get index] {
        set chain [[atomselect $id "index $indAt"] get chain]
        if {[lsearch ${l_chain} $chain] == -1} {lappend l_chain $chain}
        }
      puts $loSt "id: $id ($pdbId): chains found: ${l_chain}"
    } else {
      puts $loSt "id: $id ($pdbId): specified chains: ${l_chain}"
      }
# processing each chain
    set chainGrid 0
    foreach chain ${l_chain} {
      set l_resName $l_resNameUsr
      if {(${l_resName} == "all") || (${l_resName} == {})} {
# look for available ligand resNames
        set l_resName {}
        if {$chain == ""} {continue}
        set ligSel [atomselect $id "not (protein or nucleic) and chain $chain"]
        foreach indAt [$ligSel get index] {
          set resName [[atomselect $id "index $indAt"] get resname]
          if {[lsearch ${l_resName} $resName] == -1} {
            if {[lsearch ${exclLigName} $resName] == -1} {
              lappend l_resName $resName}}
          }
        $ligSel delete
        puts $loSt "id: $id ($pdbId): chain: $chain: ligand residue names found: ${l_resName}"
      } else {
        puts $loSt "id: $id ($pdbId): chain: $chain: specified ligand residue names: ${l_resName}"
        }
# processing each ligand resname
      foreach resName ${l_resName} {
        set l_resId $l_resIdUsr
        if {(${l_resId} == "all") || (${l_resId} == {})} {
# look for available ligand resIds
          set l_resId {}
          if {$resName == ""} {continue}
          set ligSel [atomselect $id "not (protein or nucleic) and chain $chain and resname $resName"]
          foreach indAt [$ligSel get index] {
            set resId [[atomselect $id "index $indAt"] get resid]
            if {[lsearch ${l_resId} $resId] == -1} {
              lappend l_resId $resId
              }
            }
          $ligSel delete
          puts $loSt "id: $id ([molinfo $id get name]): chain: $chain: resName: $resName: ligand resIds found: ${l_resId}"
        } else {
          puts $loSt "id: $id ([molinfo $id get name]): chain: $chain: resName: $resName: specified ligand resId: ${l_resId}"
          }
# processing each ligand resids
        foreach resId ${l_resId} {
          set writeLigand 1
          set l_ligRes [list $pdbId $chain $resName $resId]
          if {$resId == ""} {continue}
          set tmpSel [atomselect $id "chain $chain and resname $resName and resid $resId"]
          if {[$tmpSel num] < $minLigNumAtm} {
            set writeLigand 0
            if $out {puts $loSt "ligand excluded by size: $pdbId $chain $resName $resId ([$tmpSel num] atoms)"}
            }
          set testIndex [lindex [$tmpSel get index] 0]
# selecting only one altloc identifier for the ligand if specified 
          set testAltloc [[atomselect $id "index $testIndex"] get altloc]
          if {$testAltloc != "{}"} {
            set tmpSel [atomselect $id "chain $chain and resname $resName and resid $resId and altloc $testAltloc"]
            if $out {
              puts $loSt " altloc ($testAltloc) considered for: $l_ligRes"
              }
            }
# excluding user-specified ligands excluded
          foreach l_ligExclude $ll_ligExclude {
            lassign $l_ligExclude pdbIdLig chainLig resNameLig resIdLig
            if {($pdbIdLig == $pdbId) || ($pdbIdLig == ".*") || ($pdbIdLig == "\".*\"")} {
              set tmpSelLig [atomselect $id "chain $chainLig and resname $resNameLig and resid $resIdLig"]
              if {[lsearch [$tmpSelLig get index] $testIndex] >= 0} {
                set writeLigand 0
                }
              $tmpSelLig delete
              }
            }
# adding user-specified residues that are part of the receptor
          foreach l_ligRec $ll_ligRec {
            lassign $l_ligRec pdbIdLig chainLig resNameLig resIdLig
            if {($pdbIdLig == $pdbId) || ($pdbIdLig == ".*") || ($pdbIdLig == "\".*\"")} {
              set tmpSelLig [atomselect $id "chain $chainLig and resname $resNameLig and resid $resIdLig"]
              if {[lsearch [$tmpSelLig get index] $testIndex] >= 0} {
# adding user-specified ligands to the receptor selection text
                set res "(chain $chain and resname $resName and resid $resId)"
                set selTxtLigRec [concat $selTxtLigRec " or " $res]
                set writeLigand 0
                }
              $tmpSelLig delete
              }
            }
          if {($writeLigand) && ([$tmpSel num] > 0)} {
# storing identified ligand residues and writing pdbIdsFile
            if $saveDir {
              set gridName "[lindex $l_ligRes 0][lindex $l_ligRes 1]-[lindex $l_ligRes 2][lindex $l_ligRes 3]"
              set ligName "[lindex $l_ligRes 2][lindex $l_ligRes 3]-[lindex $l_ligRes 0][lindex $l_ligRes 1]"
              $tmpSel writepdb "${workPath}lig/pdb/refs/$ligName.pdb"
              $tmpSel writepdb "${workPath}lig/nopt/noh/pdb/$ligName.pdb"
              exec mkdir -p "${workPath}grid/$gridName"
              exec mkdir -p "${workPath}dock/$gridName"
# storing gridcenter.txt file with coordinates for the script gengpf.sh
              lassign [measure center $tmpSel] cx cy cz
              exec echo "[format "%.2f,%.2f,%.2f" $cx $cy $cz]" > "${workPath}grid/$gridName/gridcenter.txt"
              exec echo "$pdbId.pdbqt" > "${workPath}grid/$gridName/recname.txt"
              }
            if $out {
              puts $loSt "ligand written: ${l_ligRes} ([$tmpSel num] atoms)"}
            puts $outIds ${l_ligRes}
            set chainGrid 1
            }
          $tmpSel delete
          }
        }
# write a grid for the specific chain if no ligand was written
      if {!$chainGrid} {
        set gridName "${pdbId}${chain}"
        exec mkdir -p "${workPath}grid/$gridName"
        exec mkdir -p "${workPath}dock/$gridName"
# storing gridcenter.txt file with coordinates for the script gengpf.sh
        lassign [measure center [atomselect $id "chain $chain"] weight mass] cx cy cz
        exec echo "[format "%.2f,%.2f,%.2f" $cx $cy $cz]" > "${workPath}grid/$gridName/gridcenter.txt"
        exec echo "$pdbId.pdbqt" > "${workPath}grid/$gridName/recname.txt"
        if $out {
          puts $loSt "grid for chain written (no ligands): $gridName"
          }
        }
      }
# excluding user-specified chains from the .pdb file
    if {$recExcl != {}} {
      set selTxtLigRec "($selTxtLigRec) and not chain $recExcl"
      }
# adding user specified ligand residues that are part of the receptor
    if $saveDir {
      [atomselect $id $selTxtLigRec] writepdb "${workPath}rec/pdb/recs/${pdbId}.pdb"
      [atomselect $id $selTxtLigRec] writepdb "${workPath}rec/pdb/${pdbId}.pdb"
      }
    if $out {puts $loSt "Receptor residues included in $pdbId: $selTxtLigRec"}
    }
  if $saveDir {
    exec echo "Use 'gengpf.sh <folders>' to generate AG4 grid parameter file (gridcenter.txt and recname.txt files required).\nUse 'run_grids.sh <folders>' to run autogrid4." > "${workPath}grid/README.txt"
    exec echo "Use 'gendpf.sh <folders>' to populate with AD4 docking parameter files from ../lig/pdbqt/*.pdbqt files.\nUse 'run_docks.sh <folders>' to run autodock4." > "${workPath}dock/README.txt"
    }
  close $outIds 
  if $out {puts $loSt "output file written: $pdbIdsFile\n$procName: Done."}
  }   ;# *** lr_pdbIdsFile ***

