# 
# Supporting Information:
#
#    Profiling the interaction of 1-phenylbenzimidazoles to cyclooxygenases
#
# reference: Gómez-Castro CZ, López-Martínez M, Hernández-Pineda J,
#            Trujillo-Ferrara JG, Padilla-Martínez II. Profiling the
#            interaction of 1-phenylbenzimidazoles to cyclooxygenases.
#            J Mol Recognit. 2019; 32:e2801. https://doi.org/10.1002/jmr.2801
# 
# lrcd.tcl v-1.1.1:
#
# Tcl script for VMD to calculate the Ligand-Receptor Contact Distance (LRCD)
# parameter that compares the interaction profile of two complexes involving
# small molecules as the ligands and macromolecules as the receptors.
#
# Several variants of the LRCD parameter are implemented:
#  LRCD   - Ligand-Receptor Contact Distance.
#  LRCND  - Ligand-Receptor Contact Normalized Distance.
#  LRCDI  - Ligand-Receptor Contact Distance Index [= 1/(1 + LRCND)].
#  LRCMD  - Ligand-Receptor Contact Manhattan Distance.
#  LRCMDI - Ligand-Receptor Contact Manhattan Distance Index [= 1/(1 + LRCMD)].
#  LRCP   - Ligand-Receptor Contact Projection.
# 
# Reference: ...
#
# This script should be sourced within the program VMD either in the TkConsole,
# in the VMD's command line or from a .vmdrc file:
#
# i.e.
# # from the VMD's TkConsole:
# source lrcd.tcl
# # from the command line:
# vmd -e lrcd.tcl
#
# Instructions:
#
#   -This library was created to be used within the TkConsole of VMD, with
#    the molecules forming ligand-receptor complexes already loaded. The
#    ligand and receptor moieties of a complex may either be part of the
#    same molecule or come from different molecules, i.e., they may be
#    provided as separated files.
#
#   -This library requires prior definition by the user of the atom selections
#    specifying at least two complexes (A and B), each involving a ligand
#    (ligA and ligB) and a receptor (recA and recB) moieties. Each of this
#    four atom selections have to specify also the Id assigned by VMD when the
#    molecule containing the atom selections was loaded.
#
#   -The atom selections and Ids are stored in a Tcl array called selInfo.
#    The array keys follow the format: selInfo(<UserLabel>,<Keyword>).
#    The user specifies a label to identify a complex (<UserLabel>) and,
#    separated by a comma, one of the following keywords: 'ligId', 'ligSelTxt',
#    'recId', 'recSelTxt', and 'title'. Some examples below:
#      set selInfo(MyComplex,ligId) 0
#      set selInfo(MyComplex,ligSelTxt) "noh"
#      set selInfo(MyComplex,recId) 1
#      set selInfo(MyComplex,recSelTxt) "protein and chain A and noh"
#      set selInfo(MyComplex,title) "First lig-rec complex declared"
#
#   -Note that the ligId and recId correspond to the Ids in VMD of already
#    loaded molecules. ligSelTxt and recSelTxt are text (enclosed by "")
#    in VMD's atom selection language. See Selection Methods in the VMD
#    user's guide.
#
#   -Once declared two or more complexes in the selInfo array, the procedures
#    lrcdMat_dct, lrcdMat_cmd, and lrcdMat_lapa may be used in sequence to
#    calculate the lrcd matrix that stores lig-rec interaction measurements.
#    Then the lrcdVec_sum procedure is used to calculate the lrcd vectors
#    for each complex with sums over each interaction type. Finally, the
#    lrcd procedure uses two lrcd vectors to calculate the LRCD parameter.
#    See lrcd_test1 procedure below for an example on these steps.
#
#   -Alternatively, the procedures lrcd_pwMat and lrcd_dlgTab perform the
#    whole sequence of steps to calculate the LRCD parameter for several
#    complexes. These procedures also show the usage of lower-level procedures.
#    The lrcd_pwMat procedure calculate all pairs in a set of complexes
#    declared in the selInfo array, and specified by a list of <UserLabel>
#    identifiers. See lrcd_test2 below for an usage example. The lrcd_dlgTab
#    procedure calculate the LRCD parameter for all the complexes in a .dlg
#    file obtained from an AutoDock (version 3 or 4) molecular docking
#    calculation. In this case, the .dlg file is loaded into VMD specifying
#    PDB format. The best binding poses from the cluster analysis performed
#    by AutoDock are loaded as frames in a single molecule. Thus, the
#    procedure runs over these frames and calculate the LRCD against a
#    common reference complex. Obviously, all poses consider the same receptor.
#    See lrcd_test3 procedure below for an usage example. Several .dlg files
#    may be processed with the proper set up of the selInfo array.
#
#   -Note that the lrcd_dlgTab does not actually use any information from the
#    .dlg file, but rather only the poses of a ligand taken from a .dlg file
#    are loaded into VMD as a PDB file with several frames. Thus, the procedure
#    is independent from the docking program used to generate the complexes 
#    whenever the resulting poses may be loaded into VMD in a single molecule.
#
#   -The procedures lrcd_test* provide examples of different ways of performing
#    the LRCD calculation:
#      -lrcd_test1: performs the whole calculation for a single
#       pair of complexes taken from two PDB structures.
#      -lrcd_test2: calculates a pairwise matrix of LRCD for a set of 16
#       COX-1/2 - ligand complexes, downloading the files
#       from the repository using the internet.
#      -lrcd_test3: calculate the LRCD for the best binding poses obtained
#       from a .dlg AutoDock molecular docking calculation taken the lowest-
#       energy pose as reference for the rest.
#      -lrcd_test4: calculate the LRCD for a list of complexes or a list
#       of docking poses (.dlg file) taking as reference LRC vectors
#       pre-calculated for a set of PDB complexes (as a data base).
#
#   -In order to run the test procedures it may be neccessary to change to
#    the directory containing the lrcd.tcl script and the subdirectories
#    test1, ..., test4 with the 'cd' command in the TkConsole of VMD.



#|-lrcd.tcl :| {condText}
#|  -Tcl-language script library for VMD to calculate the LRCD parameter
#|   _ and its variants .
#|  -authors :-Carlos Z. Gómez-Castro ;
#|  -reference :
#|    -J. Mol. Recognit. 2019; 32:e2801. https://doi.org/10.1002/jmr.2801 ;
#|  -date :-2020-09-09.Wed ;
#|  -version :-1.1.1 ;
#|  -version information :
#|    -changes in this version :
#|      -new variable argument 'minLigSize' added to proc lr_pdbIdsFile :
#|        -used to exclude small residues from being considered as ligands ;
#|      -adding default filters for small ligands in proc lr_pdbIdsFile ;
#|    -finished version ;
#|  -notes from previous versions :
#|    -changes in v111 :
#|      -the output file of proc lr_pdbIdsFile is now appended instead of
#|       _ rewritten .
#|      -the 'src' argument of proc lr_pdbIdsFile is no longer variable .
#|      -the directory tree and some readme files updated to the latest scheme
#|       _ of VS scripts for ad4 (lr_pdbIdsFile) .
#|      -the procedure lr_pdbIdsFile is separated from the main script .
#|      -it is rather sourced from the main script lrcd.tcl ;
#|    -implemented excluding receptor chains in proc lr_pdbIdsFile .
#|    -the folder tree created by proc lr_pdbIdsFile is no longer
#|     _ created by default; the 'workPath' argument is required for that .
#|    -updated output in lrcdMat_dct (chain id missing) .
#|    -case-insensitive variable arguments in lr_pdbIdsFile proc .
#|    -lr_pdbIdsFile: writting 'grid' and 'dock' folder tree .
#|    -each grid folder written will contain its own 'gridcenter.txt' file .
#|    -readme files added to folder tree in ld_pdbIdsFile .
#|    -improving ligand selection/exclusion scheme in lr_pdbIdsFile :
#|      -now regular expressions may be used to specify ligands .
#|    -implemented procedure lr_pdbIdsFile :
#|      -automating detection of LR complexes in PDB structures .
#|      -cleaning structures for docking calculations .
#|      -updated lrcdVec_writeDB ;
#|    -for lrcd_dlgTab, the lrcd* values are stored into the user field .
#|    -the id of the ligand is included in the table reported by lrcd_dlgTab .
#|    -few bugs in the output of lrcd_dbTab were corrected .
#|    -in progress improving the naming of log files .
#|    -implemented and tested procedure lrcd_test4 ... .
#|    -implemented and tested procedure lrcdVec_writeDB and lrcdVec_readDB .
#|    -implemented and tested procedure lrcd_dbTab .
#|    -procedure changeDir renamed to changeVecDir .
#|    -procedures lrcdi and lrcmdi implemented .
#|    -procedures lrcd_test1, lrcd_pwMat, and lrcd_dlgTab updated .
#|    -procedure lrcmd implemented and tested :
#|      -the results seem proportional to the Cartesian distance of the
#|       _ lrcdVec of the two complexes but with higher (2-fold) values .
#|      -the normalized version (normalized lrcdVecs) is equivalent to the
#|       _ normalized lrcd .
#|      -the metric 1/(1 + lrcmd) provide greater dispersion in the
#|       _ [0,1] range, thus seems to be more convenient than the
#|       _ equivalent with the Cartesian distance ;
#|    -procedure lrcp was implemented and tested :
#|      -the values obtained follow the same trend as the normalized
#|       _ distances between lrcdVecs but are in every case very close
#|       _ to 1, thus may be useful only if the vectors are orthogonal
#|       _ (i.e. when the nature of the interactions differs greately) ;
#|    -procedure ligRecVecNorm renamed to lrcdVec_len .
#|    -procedures lrcd_pwMat and lrcd_dlg were modified to include an
#|     _ argument to specify the variant of the LRCD parameter to calculate ;
#|  -list of procedures for testing the library :
#|    -lrcd_test1 .
#|    -lrcd_test2 .
#|    -lrcd_test3 .
#|    -lrcd_test4 ;
#|  -list of procedures for calculating LR contacts (lrcdMat) :
#|    -lrcdMat_dct .
#|    -lrcdMat_cmd .
#|    -lrcdMat_lapa ;
#|  -list of procedures for calculating LR contacts vectors (lrcdVec) :
#|    -lrcdVec_sum ;
#|  -list of procedures for storing lrcdVec data sets :
#|    -lrcdVec_writeDB .
#|    -lrcdVec_readDB ;
#|  -list of procedures for calculating LRCD parameters :
#|    -lrcd .
#|    -lrcnd .
#|    -lrcdi .
#|    -lrcmd .
#|    -lrcmdi .
#|    -lrcp ;
#|  -list of procedures for calculating LRCD parameter on sets of molecules :
#|    -lrcd_pwMat .
#|    -lrcd_dlgTab .
#|    -lrcd_dbTab ;
#|  -list of miscellaneous procedures :
#|    -changeVecDir .
#|    -lrcdVec_len .
#|    -lrcdVec_norm ;
#|  -list of external procedures :
#|    -lr_pdbIdsFile.tcl ;

# public version (started from library anMol-v.0.1.4)
global lrcd_version selInfo   ;# global variables
set lrcd_version 1.1.1

# printing info when sourcing the library

puts "\n*********************************************************************"
puts "*****   Ligand-Receptor Contact Distance (LRCD) calculation   *******"
puts "***********************  version: $lrcd_version  ****************************"
puts "\n Script for calculating the LRCD parameter to compare patterns of"
puts "  ligand-receptor interactions between pairs of complexes."
puts "\n Reference: ..."
puts "   If you find the results from this script useful, please consider citing.\n"
puts "\n Quick instructions:"
puts "   1. Source the lrcd.tcl script into VMD."
puts "   2. Set up the Tcl selInfo array with ligands and receptors selections."
puts "   3. Calculate the LRCD using one of the test procedures as example."
puts "   4. See further details in the comments of the lrcd.tcl file."
puts "\n The 'lrcd_test1', 'lrcd_test2', and 'lrcd_test3' procedures may used"
puts "  to test the script and learn how to set up the calculations."
puts "  These require the directories test1, test2, and test3"
puts "  and files within to run.\n"

puts " Reading external procedures:"
if {[file exist "lr_pdbIdsFile.tcl"]} {
  source lr_pdbIdsFile.tcl
  puts "  lr_pdbIdsFile.tcl script loaded."
} elseif {[info exists tclScriptPath]} {
  if {[file exist ${tclScriptPath}lrcd-vmd.tcl/tools/lr_pdbIdsFile.tcl]} {
    source ${tclScriptPath}lrcd-vmd.tcl/tools/lr_pdbIdsFile.tcl
    puts "${tclScriptPath}lrcd-vmd.tcl/tools/lr_pdbIdsFile.tcl script sourced."
  } elseif {[file exist ${tclScriptPath}lrcd/tools/lr_pdbIdsFile.tcl]} {
    source ${tclScriptPath}lrcd/tools/lr_pdbIdsFile.tcl
    puts "${tclScriptPath}lrcd/tools/lr_pdbIdsFile.tcl script sourced."
  } else {puts "  lr_pdbIdsFile.tcl script not found."}
} else {puts "  lr_pdbIdsFile.tcl script not found."}


#|-proc lrcd_test1 {} :
#|  -performs the LRCD parameter calculation for a couple of complexes .
#|  -this procedure requires the files 1cqe.pdb and 3phg.pdb downloaded from
#|   _ the Protein Data Bank (https://www.rcsb.org/) in a directory named
#|   _ "test1" .
#|  -notes :
#|    -model 1cqe.pdb contains a homodimer of COX-1 co-crystallized with
#|     _ flurbiprofen (Flp) .
#|    -model 3pgh.pdb contains two homodimers of COX-2 co-crystallized with ;;
#|     _ flurbiprofen (Flp) ;;
proc lrcd_test1 {} {
  global selInfo lrcd_version
  puts "\nlrcd.tcl-v.${lrcd_version}: LRCD calculation of a pair of complexes."

# user specification of the lig and rec atom selections of the two complexes
# loading structure 1cqe.pdb (complexA) and setting up selInfo array
  set idA [mol new "test1/1cqe.pdb" waitfor all]
  set selInfo(complexA,ligId) $idA
  set selInfo(complexA,ligSelTxt) "resname FLP and chain A and noh"
  set selInfo(complexA,recId) $idA
  set selInfo(complexA,recSelTxt) "protein and chain A and noh"
  set selInfo(complexA,title) "Complex COX-1 - Flurbiprofen (Flp), chain A, 1CQE"

# loading structure 3pgh.pdb (complexB) and setting up selInfo array
  set idB [mol new "test1/3pgh.pdb" waitfor all]
  set selInfo(complexB,ligId) $idB
  set selInfo(complexB,ligSelTxt) "resname FLP and chain A and noh"
  set selInfo(complexB,recId) $idB
  set selInfo(complexB,recSelTxt) "protein and chain A and noh"
  set selInfo(complexB,title) "Complex COX-2 - Flurbiprofen (Flp), chain A, 3PGH"

# lrcd calculation
# calculating lrcdMat array with all 6 components for complexA
  set lrcdMat {}
  set lrcdMat [lrcdMat_dct $lrcdMat complexA 4.5]
  set lrcdMat [lrcdMat_cmd $lrcdMat complexA 4.5]
  set lrcdMat [lrcdMat_lapa $lrcdMat complexA 4.5]
# calcualting lrcdVec array for complexA
  set lrcdVecA [lrcdVec_sum $lrcdMat $selInfo(complexA,title)]

# calculating lrcdMat array with all 6 components for complexB
  set lrcdMat {}
  set lrcdMat [lrcdMat_dct $lrcdMat complexB 4.5]
  set lrcdMat [lrcdMat_cmd $lrcdMat complexB 4.5]
  set lrcdMat [lrcdMat_lapa $lrcdMat complexB 4.5]
# calcualting lrcdVec array for complexB
  set lrcdVecB [lrcdVec_sum $lrcdMat $selInfo(complexB,title)]

  puts "\nCalculatig different Lig-Rec Contact Distances between Flp-COX-1 and Flp-COX-2:"
# LRCD_(complexA,complexB)
  set lrcdRes [lrcd $lrcdVecA $lrcdVecB "LRCD_(Flp-COX-1,Flp-COX-2)" "none"]
  puts "\nLig-Rec Contact Distance (LRCD) = [format "%.1f" $lrcdRes]"
# LRCND_(complexA,complexB)
  set lrcdRes [lrcnd $lrcdVecA $lrcdVecB "LRCND_(Flp-COX-1,Flp-COX-2)" "none"]
  puts "\nLig-Rec Contact Normalized Distance (LRCND) = [format "%.3f" $lrcdRes]"
# LRCDI_(complexA,complexB)
  set lrcdRes [lrcdi $lrcdVecA $lrcdVecB "LRCDI_(Flp-COX-1,Flp-COX-2)" "none"]
  puts "\nLig-Rec Contact Distance Index (LRCDI) = [format "%.3f" $lrcdRes]"
# LRCMD_(complexA,complexB)
  set lrcdRes [lrcmd $lrcdVecA $lrcdVecB "LRCMD_(Flp-COX-1,Flp-COX-2)" "none"]
  puts "\nLig-Rec Contact Manhattan Distance (LRCMD) = [format "%.3f" $lrcdRes]"
# LRCMDI_(complexA,complexB)
  set lrcdRes [lrcmdi $lrcdVecA $lrcdVecB "LRCMDI_(Flp-COX-1,Flp-COX-2)" "none"]
  puts "\nLig-Rec Contact Manhattan Distance Index (LRCMDI) = [format "%.3f" $lrcdRes]"
# LRCP_(complexA,complexB)
  set lrcdRes [lrcp $lrcdVecA $lrcdVecB "LRCP_(Flp-COX-1,Flp-COX-2)" "none"]
  puts "\nLig-Rec Contact Projection (LRCP) = [format "%.3f" $lrcdRes]"

  puts "\nlrcd_test1: Done."
  }   ;# lrcd_test1



#|-proc lrcd_test2 {} :
#|  -performs the LRCD parameter calculation for a set of complexes using the
#|   _lrcd_pwMat procedure .
#|  -all pdb structures are attempted to be downloaded from the PDB repository ;
proc lrcd_test2 {} {
  global selInfo lrcd_version
  puts "\nlrcd.tcl-v.${lrcd_version}: Pairwise LRCD matrix of a set of complexes."
  if {![file exist test2]} {exec mkdir test2}

# user specification of the lig and rec atom selections of all the complexes.
# creating a list of the PDB Id structures to be processed.
  set pdbid [list 1diy 1eqg 1cqe 1pge 1pgf 1pgg 1fe2 1cvu 1cx2 3pgh 1pxx 4cox 6cox 1ddx]
# creating a list for the corresponding PDB residue names of the ligands
  set ligResName [list ACD IBP FLP ISF IMM IMM LAX ACD S58 FLP DIF IMN S58 PGX]

# downloading molecules and setting up selInfo array
  puts "Downloading PDB models: $pdbid"
  foreach pdb $pdbid lig $ligResName {
    set id [mol pdbload $pdb]
    set selInfo($pdb,ligId) $id
    set selInfo($pdb,ligSelTxt) "resname $lig and chain A and noh"
    set selInfo($pdb,recId) $id
    set selInfo($pdb,recSelTxt) "protein and chain A and noh"
    set selInfo($pdb,title) "Complex COX - $lig from structure $pdb"
    }

# calculating pirwise LRCD matrix for the set of complexes
  lrcd_pwMat $pdbid 4.5 "d" "test2/lrcd_pwMat"
  puts "\nlrcd_test2: Done."
  }   ;# lrcd_test2



#|-proc lrcd_test3 {} :
#|  -perform the LRCD calculation for the binding poses in an AutoDock .dlg
#|   _ file .
#|  -this procedure requires the files cox13.pdbqt and cox13-indom.dlg in
#|   _ a directory named 'test3' to work ;
proc lrcd_test3 {} {
  global selInfo lrcd_version
  puts "\nlrcd.tcl-v.${lrcd_version}: LRCD calculation for docking poses."

# user specification of the lig and rec atom selections of the complex.
# loading the receptor structure cox13.pdbqt
  set idR [mol new "test3/cox13.pdbqt" type pdb waitfor all]
  set selInfo(COX1-indom,recId) $idR
  set selInfo(COX1-indom,recSelTxt) "noh"
# NOTE: Since the cox13.pdbqt structure was preprocessed, the VMD macros for
#       atom selection such as 'protein' will not work. In this case all
#       receptor atoms except for hydrogens are considered.

# loading the .dlg file with the calculated ligand poses
  set idL [mol new "test3/cox13-indom.dlg" type pdb waitfor all]
  set selInfo(COX1-indom,ligId) $idL
  set selInfo(COX1-indom,ligSelTxt) "noh"
  set selInfo(COX1-indom,title) "Complex COX-1 - Indomethacin (Imn), docking poses."

# calculating the LRCDMDI parameter for all the ligand poses
  mol top $idL
  animate goto start   ;# move to the firts frame to use it as reference
  lrcd_dlgTab [list COX1-indom] COX1-indom 4.5 "mdi" "test3/lrcmdi_dlgTab"
  puts "\nlrcd_test3: Done."
  }   ;# lrcd_test3



#|-lrcd_test4 {} :
#|  -performs the LRCD parameter calculation of set of docking poses against
#|   _ a data base of complexes using pre-calculated LRC vectors ;
proc lrcd_test4 {} {
  global selInfo lrcd_version
  puts "\nlrcd.tcl-v.${lrcd_version}: LRCD for a pre-calculated ref data base."

# user specification of the docking complex to be calculated (cox13-indom)
# loading the receptor structure cox13.pdbqt
  set idR [mol new "test4/cox13.pdbqt" type pdb waitfor all]
  set selInfo(COX1-indom,recId) $idR
  set selInfo(COX1-indom,recSelTxt) "noh"
# NOTE: Since the cox13.pdbqt structure was preprocessed, the VMD macros for
#       atom selection such as 'protein' will not work. In this case all
#       receptor atoms except for hydrogens are considered.

# loading the .dlg file with the calculated ligand poses
  set idL [mol new "test4/cox13-indom.dlg" type pdb waitfor all]
  set selInfo(COX1-indom,ligId) $idL
  set selInfo(COX1-indom,ligSelTxt) "noh"
  set selInfo(COX1-indom,title) "Complex COX-1 - Indomethacin (Imn), docking poses."

# loading the pre-calculated data set of lrcdVec of PDB reference complexes
# loading data for 24 COX-1 complexes, and then adding 50 COX-2 complexes
  set lrcdVecArr {}   ;# initilizing
  set lrcdVecArr [lrcdVec_readDB "test4/lrcdVecDB_cox1-4.5.txt" $lrcdVecArr]
  set lrcdVecArr [lrcdVec_readDB "test4/lrcdVecDB_cox2-4.5.txt" $lrcdVecArr]

# calculating the LRCDI parameter for all the ligand poses
  lrcd_dbTab [list COX1-indom] $lrcdVecArr 4.5 "d" "test4/lrcmdi_dbTab"
  puts "\nlrcd_test4: Done."
  }   ;# lrcd_test4




#|-proc lrcdMat_dct {lrcdMat cmplxSelId cutoff {loSt stdout}} :
#|  -incorporates the distances by contact type (DCT) into the lrcdMat array .
#|  -adds the keys 'pp', 'pn', 'np' and 'nn' to lrcdMat .
#|  -each key contains a list of lists with format {{indAtL indAtR dist} ...} .
#|  -return a list formatted to create an array with the matrix information :
#|    -array format :
#|      -arrayname(key) {{index index distance} ...} ;
#|  -definition of the array lrcdMat ... .
#|  -arguments :
#|    -lrcdMat :-ligand-receptor contact distance (lrcd) matrix .
#|      -acceptable values :
#|        -array with a lrcdMat in tcl list format :
#|          -if the array already have DCT information the old distances will
#|           _ be overwritten ;
#|        -an empty list ({}) or ("") ;;
#|    -cmplxSelId :-selId for the selInfo array with lig-rec complex info .
#|      -selInfo keys required (selInfo(cmplxSelId,key)) :
#|        -title, ligId, ligSelTxt, recId and recSelTxt .
#|        -for the selTxt 'none' is acceptable to override the calculation ;;
#|    -cutoff :-cutoff lig-rec distance limit to be considered in the lrcdMat .
#|      -acceptable values :
#|        -positive real values in Angstroms .
#|        -a value of 4.5 yields good results ;;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none' :-will avoid any output ;;;;
#|  -notes :
#|    -requires access to the selInfo global array .
#|    -this procedure replaces ligRecDistMat .
#|    -a showSelInfo procedure for lig-rec complex info would be useful ;;
proc lrcdMat_dct {lrcdMat cmplxSelId cutoff {loSt stdout}} {
# extract and check input info
  global selInfo
  array set lrcdArr $lrcdMat
  if {$loSt == "none"} {set out 0} else {set out 1}
  set ligId $selInfo($cmplxSelId,ligId)
  set ligTxt $selInfo($cmplxSelId,ligSelTxt)
  set recId $selInfo($cmplxSelId,recId)
  set recTxt $selInfo($cmplxSelId,recSelTxt)
  if {$out} {
    puts $loSt "\n+++ ligand-receptor contact distance matrix +++"
    puts $loSt "+++ distances by contact type (types: pp pn np nn) +++"
    puts $loSt "\ncutoff: $cutoff"
    puts $loSt "selInfo selId: $cmplxSelId   title: $selInfo($cmplxSelId,title)"
    puts $loSt "ligand: id: $ligId; sel: $ligTxt"
    puts $loSt "receptor: id: $recId; sel: $recTxt"
    }
  if {($ligTxt == "none") || ($recTxt == "none")} {return $lrcdMat}
# calculate lig-rec distances by contact type within cutoff and feed lrcdMat
  set lig [atomselect $ligId $ligTxt]
  set rec [atomselect $recId $recTxt]
  foreach ct {pp pn np nn} {set lrcdArr($ct) {}}
  foreach lAtInd [$lig get index] {
    foreach rAtInd [$rec get index] {
      set atl [list [list $lAtInd $ligId] [list $rAtInd $recId]]
      set dist [measure bond $atl]
      if {$dist <= $cutoff} {
        set lAtS [atomselect $ligId "index $lAtInd"]
        set rAtS [atomselect $recId "index $rAtInd"]
        set lAtType [string index [$lAtS get name] 0]
        set rAtType [string index [$rAtS get name] 0]
        if {($lAtType == "O")||($lAtType == "N")||($lAtType == "S")\
            ||($lAtType == "P")} {   ;# ligand atom polar
          if {($rAtType == "O")||($rAtType == "N")||($rAtType == "S")\
              ||($rAtType == "P")} {   ;# receptor atom also polar
            lappend lrcdArr(pp) [list $lAtInd $rAtInd $dist]
          } else {   ;#   receptor atom non-polar
            lappend lrcdArr(pn) [list $lAtInd $rAtInd $dist]
            }
        } else {   ;# lig atom non-polar
          if {($rAtType == "O")||($rAtType == "N")||($rAtType == "S")\
              ||($rAtType == "P")} {   ;# receptor atom polar
            lappend lrcdArr(np) [list $lAtInd $rAtInd $dist]
          } else {   ;# receptor atom also non-polar
            lappend lrcdArr(nn) [list $lAtInd $rAtInd $dist]
            }
          }
        $lAtS delete
        $rAtS delete
        }
      }
    }
  $lig delete
  $rec delete
# output results summary
  if {$out} {
    foreach cType [list pp pn np nn] {
      set nCont [llength $lrcdArr($cType)]
      switch $cType {
        pp {puts $loSt "\n $nCont polar - polar L-R contacts:"}
        pn {puts $loSt "\n $nCont polar - non-polar L-R contacts:"}
        np {puts $loSt "\n $nCont non-polar - polar L-R contacts:"}
        nn {puts $loSt "\n $nCont non-polar - non-polar L-R contacts:"}
        }
      puts $loSt " ligand atom \t receptor atom \t distance (Angstroms)"
      set lrcdArr($cType) [lsort -index 2 -real $lrcdArr($cType)]
      foreach iid $lrcdArr($cType) {
        lassign $iid indL indR dist
        set lAtS [atomselect $ligId "index $indL"]
        set rAtS [atomselect $recId "index $indR"]
        puts $loSt "[molinfo $ligId get name].[$lAtS get resname][$lAtS get resid].[$lAtS get name].[$lAtS get index]\t[molinfo $recId get name].[$rAtS get chain].[$rAtS get resname][$rAtS get resid].[$rAtS get name].[$rAtS get index]\t$dist"
        $lAtS delete
        $rAtS delete
        }
      }
    puts $loSt "\nlrcdMat_dct: Done."
    }
  return [array get lrcdArr]
  }   ;# *** lrcdMat_dct ***

#|-proc lrcdMat_cmd {lrcdMat cmplxSelId cutoff {loSt stdout}} :
#|  -incorporates the distance between centers of mass (COM) from a ligand and
#|   _ a receptor into the lrcdMat array .
#|  -the COM-distance is scaled by the cutoff distance .
#|  -return a list formatted to create an array with the matrix information :
#|    -array format :
#|      -lrcdMat(key) {"COM-lig" "COM-rec" distance} ;
#|    -keys (names) and their content :
#|      -pp :-lists of indices and distance between polar atoms in both the
#|       _ ligand and in the receptor within the cutoff distance;
#|      -pn :-lists of indices and distance between polar atoms in the ligand
#|       _ and non-polar atoms in the receptor within the cutoff distance ;
#|      -np :-lists of indices and distance between non-polar atoms in the
#|       _ ligand and polar atoms in the receptor within the cutoff distance ;
#|      -nn :-lists of indices and distance between non-polar atoms in both
#|       _ ligand and receptor atoms within the cutoff distance ;
#|      -cmd :-distance between COM's from ligand and receptor .
#|        -the indices refer to the first atom on each L or R selections ;;
#|    -the output list can be converted into an array with the command
#|     _ 'array set arrVar $outList' ;
#|  -arguments :
#|    -lrcdMat :-lig-rec contact distance matrix (array) in tcl list fotmat ;
#|    -cmplxSelId :-selId for the selInfo array with lig-rec complex info .
#|      -selInfo keys required (selInfo(cmplxSelId,key)) :
#|        -title :-user information or description of the complex ;
#|        -ligId :-Id of the molecule treated as ligand ;
#|        -ligSelTxt :-atom selection text refering to the ligand molecule ;
#|        -recId  :-Id of the molecule treated as receptor ;
#|        -recSelTxt :-atom selection text refering to the receptor moiety ;;
#|      -if ligSelTxt or refSelTxt is set to 'none' no action is performed
#|       _ and the original lrcd matrix is returned ;
#|    -cutoff :-maximum distance limit for L-R contacts to be considered .
#|      -this value is also used to scale the distance between COM's ;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none':-will avoid any output ;;;;
#|  -notes :
#|    -the procedure requires acces to the global selInfo array ;;
proc lrcdMat_cmd {lrcdMat cmplxSelId cutoff {loSt stdout}} {
# extract and check input info
  global selInfo
  array set lrcdArr $lrcdMat
  if {$loSt == "none"} {set out 0} else {set out 1}
  set ligId $selInfo($cmplxSelId,ligId)
  set ligTxt $selInfo($cmplxSelId,ligSelTxt)
  set recId $selInfo($cmplxSelId,recId)
  set recTxt $selInfo($cmplxSelId,recSelTxt)
  if {$out} {
    puts $loSt "\n+++ ligand-receptor contact distance matrix +++"
    puts $loSt "+++ distance between center of mass (COM) of lig and rec +++"
    puts $loSt "\ncutoff: $cutoff"
    puts $loSt "selInfo selId: $cmplxSelId   title: $selInfo($cmplxSelId,title)"
    puts $loSt "ligand: id: $ligId; sel: $ligTxt"
    puts $loSt "receptor: id: $recId; sel: $recTxt"
    }
  if {($ligTxt == "none") || ($recTxt == "none")} {return $lrcdMat}
# calculate distance between COMs
  set ligSel [atomselect $ligId $ligTxt]
  lassign [measure center $ligSel weight mass] xl yl zl
  set recSel [atomselect $recId $recTxt]
  lassign [measure center $recSel weight mass] xr yr zr
  set comD [expr {sqrt(($xl - $xr)**2 + ($yl - $yr)**2 + ($zl - $zr)**2)}]
  $ligSel delete
  $recSel delete
# report results
  if {$out} {
    puts $loSt "\nCOM-lig: ($xl $yl $zl)"
    puts $loSt "COM-rec: ($xr $yr $zr)"
    puts $loSt "distance between COMs (A): $comD; scaled: [expr $comD*$cutoff]"
    puts $loSt "\nlrcdMat_cmd: Done."
    }
  set lrcdArr(cmd) [list [list "COM-lig" "COM-rec" [expr $comD*$cutoff]]]
  return [array get lrcdArr]
  } ;#   *** lrcdMat_cmd ***

#|-proc lrcdMat_lapa {lrcdMat cmplxSelId cutoff {loSt stdout}} :
#|  -incorporates the length of the arc between the principal axes of inertia
#|   _ of a ligand and a receptor into an array lrcdMat .
#|  -for each principal axis of inertia, the angle between the lig and rec axes 
#|   _ is considered to calculate the length of the circular arc with radius
#|   _ equal to the cutoff distance .
#|  -a direction for each principal axis is conventionally choosen as the
#|   _ "most massive" sides for both lig and rec molecules from their COM .
#|  -three arc lengths are reported as a measure of the relative
#|   _ orientations between the ligand and the receptor .
#|  -returns a list formatted to create an array with the matrix information :
#|    -array format :
#|      -arrayname(key) {{PAIi-lig PAIi-rec length} ...} ;
#|    -keys (names) and their content :
#|      -pp :-lists of indices and distance between polar atoms in both the
#|       _ ligand and in the receptor within the cutoff distance;
#|      -pn :-lists of indices and distance between polar atoms in the ligand
#|       _ and non-polar atoms in the receptor within the cutoff distance ;
#|      -np :-lists of indices and distance between non-polar atoms in the
#|       _ ligand and polar atoms in the receptor within the cutoff distance ;
#|      -nn :-lists of indices and distance between non-polar atoms in both
#|       _ ligand and receptor atoms within the cutoff distance ;
#|      -cmd :-squared distance between COM's from ligand and receptor .
#|        -dummy indices are included to have conserve the format .
#|        -the indices refer to the first atom on each L or R selections ;
#|      -arc :-lengths of the arcs formed between
#|       _ the principal axes of inertia of the ligand and the receptor .
#|        -the radius of the circle is the same as the cutoff argument ;;
#|    -the output list can be converted into an array with the command
#|     _ 'array set arrVar $outList' ;
#|  -arguments :
#|    -lrcdMat :-lig-rec contact distance matrix (array) in tcl list fotmat ;
#|    -cmplxSelId :-selId for the selInfo array with lig-rec complex info .
#|      -selInfo keys required (selInfo(cmplxSelId,key)) :
#|        -title :-user information or description of the complex ;
#|        -ligId :-Id of the molecule treated as ligand ;
#|        -ligSelTxt :-atom selection text refering to the ligand molecule ;
#|        -recId  :-Id of the molecule treated as receptor ;
#|        -recSelTxt :-atom selection text refering to the receptor moiety ;
#|      -if ligSelTxt or refSelTxt is set to 'none' :
#|        -no action is performed and the original lrcd matrix is returned ;;;
#|    -cutoff :-maximum distance limit for L-R contacts to be considered .
#|      -this value is also used to scale the length of the arcs between
#|       _ principal axes of inertia ;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file .
#|        -'none':-will avoid any output ;;;;
#|  -notes :
#|    -by symmetry, different orientations of a molecule may yield similar
#|     _ lengths of the arc between the angles of one of its principal axes
#|     _ with respect to the principal axes from other molecule .
#|    -the procedure requires access to the global selInfo array ;;
proc lrcdMat_lapa {lrcdMat cmplxSelId cutoff {loSt stdout}} {
# extract and check input info
  global selInfo
  array set lrcdArr $lrcdMat
  if {$loSt == "none"} {set out 0} else {set out 1}
  set ligId $selInfo($cmplxSelId,ligId)
  set ligTxt $selInfo($cmplxSelId,ligSelTxt)
  set recId $selInfo($cmplxSelId,recId)
  set recTxt $selInfo($cmplxSelId,recSelTxt)
  if {$out} {
    puts $loSt "\n+++ ligand-receptor contact distance matrix +++"
    puts $loSt "+++ length of the arc between principal axes of lig and rec +++"
    puts $loSt "\ncutoff: $cutoff"
    puts $loSt "selInfo selId: $cmplxSelId   title: $selInfo($cmplxSelId,title)"
    puts $loSt "ligand: id: $ligId; sel: $ligTxt"
    puts $loSt "receptor: id: $recId; sel: $recTxt\n"
    }
  if {($ligTxt == "none") || ($recTxt == "none")} {return $lrcdMat}
# calculate length of the arc between each of the principal axes of inertia
  set lrcdArr(arc) {}
  set ligSel [atomselect $ligId $ligTxt]
  set recSel [atomselect $recId $recTxt]
# principal axis of inertia
  lassign [measure inertia $ligSel] comL paL
  lassign [measure inertia $recSel] comR paR
  foreach pAxisIL $paL pAxisIR $paR axis {1 2 3} {
# find the most massive direction on each axis
    set pAxisIL [changeVecDir $pAxisIL $comL $ligId $ligTxt]
    set pAxisIR [changeVecDir $pAxisIR $comR $recId $recTxt]
    set angleRad [expr {acos([vecdot $pAxisIL $pAxisIR])}]
    set arcL [expr $angleRad*$cutoff]
    lappend lrcdArr(arc) [list "PAI$axis-lig" "PAI$axis-rec" $arcL]
    if {$out} {puts $loSt "Angle between principal axis $axis L-R: [expr {$angleRad*180.0/3.141592}] deg; length of arc: $arcL A"}
    }
  $ligSel delete
  $recSel delete
  if {$out} {puts $loSt "\nlrcdMat_lapa: Done."}
  return [array get lrcdArr]
  } ;# *** lrcdMat_lapa ***

# used by the lrcdMat_lapa procedure
proc changeVecDir {vec com id selTxt} {
  lassign $vec vx vy vz
  set d [vecdot $vec $com]
  set totMass [eval vecadd [[atomselect $id $selTxt] get mass]]
  set posDirMass [eval vecadd [[atomselect $id "$selTxt and ($vx*x + $vy*y + $vz*z > $d)"] get mass]]
  if {$posDirMass < [expr {$totMass - $posDirMass}]} {
    return [vecscale -1.0 $vec]} else {return $vec}
  }

#|-proc lrcdVec_sum {lrcdMat {infoM ""} {loSt stdout}} :
#|  -calculates a lrcd vector where each component contain the sum of all the
#|   _ measurements contained in each category (key) of the lrcdMat .
#|  -returns a list formated to create an array with the sums of each category .
#|  -arguments :
#|    -lrcdMat :-lig-rec contact distance matrix (array) in tcl list fotmat .
#|      -obtained from the output of one or more of the procedures lrcdMat_* ;
#|    -infoM :-string with information about the system .
#|      -default value :-empty string ("") ;;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none':-will avoid any output ;;;;
proc lrcdVec_sum {lrcdMat {infoM ""} {loSt stdout}} {
  array unset lrcdArr
  array set lrcdArr $lrcdMat
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {$out} {
    puts $loSt "\n+++ ligand-receptor contact distance vector +++"
    puts $loSt "+++ sum of values for each component +++"
    puts $loSt "\nuser info: $infoM"
    puts $loSt "categories in the lrcd matrix: [array names lrcdArr]"
    puts $loSt "\n type \t N \t dist sum (Angstroms)"
    }
  foreach type [array names lrcdArr] {
    set distList {}
    foreach contDist $lrcdArr($type) {
      lappend distList [lindex $contDist 2]
      }
    set lrcvsd($type) 0.0
    foreach dist $distList {
      set lrcvsd($type) [expr {$lrcvsd($type) + $dist}]
      }
    if {$out} {puts $loSt "$type\t[llength $distList]\t$lrcvsd($type)"}
    }
  if {$out} {puts $loSt "\nlrcdVec_sum: Done."}
  return [array get lrcvsd]
  }   ;# *** lrcdVec_sum ***


#|-lrcdVec_writeDB {pdbIdsFile {cutoff 4.5} {lrcdVecFile "lrcdVecDB.txt"}
#|                  _ {lrcdMatL {dct cmd lapa}} {loSt stdout}} :
#|  -calculate the LRCD vector for a list of PBD complexes and store them into
#|   _ a text file .
#|  -the proc lrcdVec_readDB may be used to recover a set of stored vectors .
#|  -the list of complexes is read from a text file with four columns
#|   _ separated by space characters .
#|  -the pdbIdsFile may be automatically generated with proc lr_pdbIdsFile .
#|  -All PDB models will be downloaded from the PDB repository .
#|  -currently, only vectors of type lrcdVec_sum are calculated, in future
#|   _ versions different types may be incorporated .
#|  -arguments :
#|    -pdbIdsFile :
#|      -name of the text file with pdbIds, chains, resNames, and resId
#|       _ identifiers (as generated by proc lr_pdbIdsFile) .
#|      -format (4 columns) :
#|        -(1) 4-letter PDBID .
#|        -(2) 1-letter chain identifier (generally A) .
#|        -(3) 3-letter PDB residue name for the ligand .
#|        -(4) 1- to 4-letter residue seq. number (resid) (new in v.1.0.6) ;
#|    -cutoff :
#|      -maximum distance limit (in Angstroms) to consider contacts and
#|       _ to scale lrcd parameters .
#|      -default value :-4.5 ;;
#|    -lrcdVecFile :
#|      -name of the text file that will be created to store the calculated
#|       _ lrcd vectors .
#|      -default value :-"lrcdVecDB.txt" ;;
#|    -lrcdMatL :
#|      -list of flags indicating which lrcdMat procedures are to be
#|       _ considered in the creation of the lrcd matrices .
#|      -accepatble values :
#|        -a tcl list with one or more of the flags 'cdt', 'cmd' and 'lapa' .
#|        -other flags may be incorporated in future versions ;
#|      -default value :-{dct cmd lapa} ;;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none':-will avoid any output ;;;;;
proc lrcdVec_writeDB {pdbIdsFile {cutoff 4.5} \
                      {lrcdVecFile "lrcdVecDB_4.5.txt"} \
                      {lrcdMatL {dct cmd lapa}} {loSt stdout}} {
  global selInfo
  set procName [lindex [info level 0] 0]
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {![file exist $pdbIdsFile]} {
    if {$out} {puts "$procName: Error: File not found: $pdbIdsFile"}
    return
    }
  if {$out} {
    puts $loSt "\n+++ $procName: Storing LRCD vectors for a set of PDB complexes +++"
    puts $loSt "\n Reading PDB list from pdbIdsFile: $pdbIdsFile"
    puts $loSt " LRCD vector DB text file to be created: $lrcdVecFile"
    puts $loSt " cutoff: $cutoff"
    }
  set outDB [open $lrcdVecFile w]
  set pdbIds [split [read [open $pdbIdsFile r]] \n]
  foreach pdbId $pdbIds {
    array unset lrcdVec
    set lrcdMat {}
    if {$pdbId == ""} {break}
    lassign $pdbId pdb chain name resid
    if {$out} {
      puts $loSt "Calculating lrcdVec: pdbId: $pdb; chain: $chain; resname: $name; resid: $resid"
      }
# downloading PDB and setting up selInfo array
    set id [mol pdbload $pdb]
    set selInfo(tmpPDB,recId) $id
    set selInfo(tmpPDB,ligId) $id
    set selInfo(tmpPDB,recSelTxt) "protein and noh"
    set selInfo(tmpPDB,ligSelTxt) "chain $chain and resname $name and resid $resid and noh"
    set selInfo(tmpPDB,title) "pdbId: $pdb; chain: $chain resname: $name; resid: $resid"
# calculating LR contacts (vector components)
    foreach flag $lrcdMatL {
      set lrcdMat [lrcdMat_$flag $lrcdMat tmpPDB $cutoff $loSt]
      }
# calculating the lrcdVec and writting to the output file
    set lrcdVec([join $pdbId "-"]) [lrcdVec_sum $lrcdMat $selInfo(tmpPDB,title) $loSt]
    puts $outDB [array get lrcdVec]
    mol delete $id
    }
  close $outDB
  if {$out} {
    puts $loSt "\nFile written: $lrcdVecFile"
    puts $loSt "\nlrcdVec_writeDB: Done."
    }
  }   ;# lrcdVec_writeDB


#|-proc lrcdVec_readDB {} :
#|  -reads from a text file (generated by the proc lrcdVec_writeDB) stored
#|   _ lrcdVec data pre-calculated for a set of PDB complexes .
#|  -returns a string to build a Tcl array that include the whole data set
#|   _ in the form lrcdVecDB(<complexId>) where each element contains a lrcd
#|   _ vector of type lrcdVec_sum by default .
#|  -no molecules are loaded or processed .
#|  -arguments :
#|    -lrcdVecFile :
#|      -name of the text file containing an array of lrcdVec arrays as
#|       _ genreated by the lrcdVec_writeDB procedure ;
#|    -lrcdVecDB :
#|      -string to generate a Tcl array with lrcdVec arrays .
#|      -serves to incorporate to a previous collection of lrcdVec arrays
#|       _ new elements read from a file .
#|      -default value :-{} ;;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none':-will avoid any output ;;;;;
proc lrcdVec_readDB {lrcdVecFile {lrcdVecDB {}} {loSt stdout}} {
  set pName [lindex [info level 0] 0]
  array set lrcdVecs $lrcdVecDB
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {![file exist $lrcdVecFile]} {
    if {$out} {puts $loSt "\n$pName: Error: File not found: $lrcdVecFile"}
    return
    }
  if {$out} {
    puts $loSt "\n+++ $pName: Read LRCD vectors for a set of PDB complexes +++"
    }
  set lrcdDB [split [read [open $lrcdVecFile r]] \n]
  set nComp 0
  foreach comp $lrcdDB {
    if {$comp == ""} {continue}
    array set lrcdVecs $comp
    incr nComp
    }
  if {$out} {
    puts $loSt "\n $nComp records were read from file: $lrcdVecFile."
    puts $loSt " total number of records: [llength [array names lrcdVecs]]"
    puts $loSt " complexIDs: [array names lrcdVecs]"
    puts $loSt "\n$pName: Done."
    }
  return [array get lrcdVecs]
  }   ;# lrcdVec_readDB


#|-lrcd {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} :
#|  -calculates the LRCD parameter between two lrcd vectors comparing the
#|   _ interaction pattern of two ligand-receptor complexes .
#|  -returns the resulting LRCD value and reports the two lrcd vectors .
#|  -in contras to proc lrcnd, the lrcd vector are not normalized here .
#|  -definitions :
#|    -lrcd :-square root of the sum of squared differences between the
#|     - components of two lig-rec contact distance vectors (lrcdVec) ;
#|    -lrcdVec :-vector (array) characterizing the contacts of
#|     _ a ligand with a receptor (complex L-R) ;;
#|  -arguments :
#|    -lrcdVec1 :-L-R contact vector for reference complex .
#|      -list formatted to build an array as returned by ligRecContVec_SD ;
#|    -lrcdVec2 :-L-R contact vector for the complex to compare to ref .
#|      -list formatted to build an array as returned by ligRecContVec_SD ;
#|    -infoV :-information string about the systems .
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none':-will avoid any output ;;;;;
proc lrcd {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} {
  array unset lrcdArr1
  array unset lrcdArr2
  array set lrcdArr1 $lrcdVec1
  array set lrcdArr2 $lrcdVec2
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {[llength [array names lrcdArr1]] != [llength [array names lrcdArr2]]} {
    if {$out} {puts $loSt "lrcd: Error: Incompatible LRCD vectors"}
    return "-"
    }
  if {$out} {
    puts $loSt "\n+++ligand-receptor contact distance (LRCD) between two complexes +++"
    puts $loSt "\nuser info: $infoV"
    puts $loSt "types considered in the L-R complex 1: [array names lrcdArr1]"
    puts $loSt "types considered in the L-R complex 2: [array names lrcdArr2]"
    puts $loSt "\n component \t complex 1 \t complex 2 \t cmplx 2 - cmplx 1"
    }
  set sum 0.0
  set n 0
  foreach type [array names lrcdArr1] {
    set dif [expr {$lrcdArr2($type) - $lrcdArr1($type)}]
    if {$out} {puts $loSt "$type\t$lrcdArr1($type)\t$lrcdArr2($type)\t$dif"}
    set sum [expr {$sum + $dif*$dif}]
    incr n
    }
  set lrcv [expr {sqrt($sum)}]
  if {$out} {
    puts $loSt "\nlig-rec contact distance: $lrcv"
    puts $loSt "\nlrcd: Done."
    }
  return $lrcv
  }   ;# *** lrcd ***



#|-lrcnd {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} :
#|  -calculates the LRCND parameter between two lrcd vectors comparing the
#|   _ interaction pattern of two ligand-receptor complexes .
#|  -the lrcd vectors are normalized before calculating the Cartesian distance .
#|  -returns the resulting LRCD value and reports the two lrcd vectors .
#|  -definitions :
#|    -lrcnd :-square root of the sum of squared differences between the
#|     - components of two lig-rec contact distance vectors (lrcdVec) ;
#|    -lrcdVec :-vector (array) characterizing the contacts of
#|     _ a ligand with a receptor (complex L-R) ;;
#|  -arguments :
#|    -lrcdVec1 :-L-R contact vector for reference complex .
#|      -list formatted to build an array as returned by ligRecContVec_SD ;
#|    -lrcdVec2 :-L-R contact vector for the complex to compare to ref .
#|      -list formatted to build an array as returned by ligRecContVec_SD ;
#|    -infoV :-information string about the systems .
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none':-will avoid any output ;;;;;
proc lrcnd {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} {
  set lrcdVec1 [lrcdVec_norm $lrcdVec1 "none"]
  set lrcdVec2 [lrcdVec_norm $lrcdVec2 "none"]
  array unset lrcdArr1
  array unset lrcdArr2
  array set lrcdArr1 $lrcdVec1
  array set lrcdArr2 $lrcdVec2
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {[llength [array names lrcdArr1]] != [llength [array names lrcdArr2]]} {
    if {$out} {puts $loSt "lrcd: Error: Incompatible LRCD vectors"}
    return "-"
    }
  if {$out} {
    puts $loSt "\n+++ligand-receptor contact normalized distance (LRCND) between two complexes +++"
    puts $loSt "\nuser info: $infoV"
    puts $loSt "types considered in the L-R complex 1: [array names lrcdArr1]"
    puts $loSt "types considered in the L-R complex 2: [array names lrcdArr2]"
    puts $loSt "\n component \t complex 1 \t complex 2 \t cmplx 2 - cmplx 1"
    }
  set sum 0.0
  foreach type [array names lrcdArr1] {
    set dif [expr {$lrcdArr2($type) - $lrcdArr1($type)}]
    if {$out} {puts $loSt "$type\t$lrcdArr1($type)\t$lrcdArr2($type)\t$dif"}
    set sum [expr {$sum + $dif*$dif}]
    }
  set lrcv [expr {sqrt($sum)}]
  if {$out} {
    puts $loSt "\nlig-rec contact normalized distance: $lrcv"
    puts $loSt "\nlrcnd: Done."
    }
  return $lrcv
  }   ;# *** lrcnd ***


#|-proc lrcdi {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} :
#|  -calculates the LRCNDI parameter between two lrcd vectors comparing the
#|   _ interaction pattern of two ligand-receptor complexes .
#|  -the lrcd vectors are normalized before calculating the Cartesian distance .
#|  -returns the resulting LRCD value and reports the two lrcd vectors .
#|  -definitions :
#|    -lrcdi :-LRCDI = 1/(1 + LRCND) ;
#|    -lrcnd :-square root of the sum of squared differences between the
#|     _ components of two normalized LR contact distance vectors (lrcdVec) ;
#|    -lrcdVec :-vector (array) characterizing the contacts of
#|     _ a ligand with a receptor (complex L-R) ;;
#|  -arguments :
#|    -lrcdVec1 :-L-R contact vector for reference complex .
#|      -list formatted to build an array as returned by ligRecContVec_SD ;
#|    -lrcdVec2 :-L-R contact vector for the complex to compare to ref .
#|      -list formatted to build an array as returned by ligRecContVec_SD ;
#|    -infoV :-information string about the systems .
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none':-will avoid any output ;;;;;
proc lrcdi {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} {
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {$out} {puts $loSt "\n+++ligand-receptor contact distance index (LRCDI) between two complexes +++"}
  set lrcdival [expr {1/(1 + [lrcnd $lrcdVec1 $lrcdVec2 $infoV $loSt])}]
  if {$out} {
    puts $loSt "\nlig-rec contact distance index: $lrcdival"
    }
  return $lrcdival
  }

#|-proc lrcmd {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} :
#|  -calculates the Manhattan-distance variant of the LRCD parameter ;
#|  -returns the resulting LRCD-MD value and reports the two lrc vectors .
#|  -the lrc vectors are allways normalized .
#|  -definitions :
#|    -lrcmd :-absolute value of the differences between the
#|     - components of two lig-rec contact distance vectors (lrcdVec) ;
#|    -lrcdVec :-vector (array) characterizing the contacts of
#|     _ a ligand with a receptor (complex L-R) ;;
#|  -arguments :
#|    -lrcdVec1 :-L-R contact vector for reference complex .
#|      -list formatted to build an array as returned by lrcdVec_sum ;
#|    -lrcdVec2 :-L-R contact vector for the complex to compare to ref .
#|      -list formatted to build an array as returned by lrcdVec_sum ;
#|    -infoV :-information string about the systems .
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none':-will avoid any output ;;;;;
proc lrcmd {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} {
  set lrcdVec1 [lrcdVec_norm $lrcdVec1 $loSt]
  set lrcdVec2 [lrcdVec_norm $lrcdVec2 $loSt]
  array unset lrcdArr1
  array unset lrcdArr2
  array set lrcdArr1 $lrcdVec1
  array set lrcdArr2 $lrcdVec2
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {[llength [array names lrcdArr1]] != [llength [array names lrcdArr2]]} {
    if {$out} {puts $loSt "lrcmd: Error: Incompatible LRCD vectors"}
    return "-"
    }
  if {$out} {
    puts $loSt "\n+++ligand-receptor contact Manhattan distance (LRCMD) between two complexes +++"
    puts $loSt "\nuser info: $infoV"
    puts $loSt "types considered in the L-R complex 1: [array names lrcdArr1]"
    puts $loSt "types considered in the L-R complex 2: [array names lrcdArr2]"
    puts $loSt "\n component \t complex 1 \t complex 2 \t cmplx 2 - cmplx 1"
    }
  set sum 0.0
  foreach type [array names lrcdArr1] {
    set dif [expr {$lrcdArr2($type) - $lrcdArr1($type)}]
    if {$out} {puts $loSt "$type\t$lrcdArr1($type)\t$lrcdArr2($type)\t$dif"}
    set sum [expr {$sum + abs($dif)}]
    }
  if {$out} {
    puts $loSt "\nlig-rec contacts (Mannhatan) distance: $sum"
    puts $loSt "\nlrcmd: Done."
    }
  return $sum
  }   ;# lrcmd


#|-proc lrcmdi {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} :
#|  -calculates the Manhattan-distance index variant of the LRCD parameter ;
#|  -returns the resulting LRCD-MD value and reports the two lrc vectors .
#|  -the lrc vectors are allways normalized .
#|  -definitions :
#|    -lrcmdi :-LRCMDI = 1/(1 + LRCMD) ;
#|    -lrcmd :-absolute value of the differences between the
#|     - components of two lig-rec contact distance vectors (lrcdVec) ;
#|    -lrcdVec :-vector (array) characterizing the contacts of
#|     _ a ligand with a receptor (complex L-R) ;;
#|  -arguments :
#|    -lrcdVec1 :-L-R contact vector for reference complex .
#|      -list formatted to build an array as returned by lrcdVec_sum ;
#|    -lrcdVec2 :-L-R contact vector for the complex to compare to ref .
#|      -list formatted to build an array as returned by lrcdVec_sum ;
#|    -infoV :-information string about the systems .
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none':-will avoid any output ;;;;;
proc lrcmdi {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} {
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {$out} {puts $loSt "\n+++ligand-receptor contact Manhattan distance index (LRCMDI) between two complexes +++"}
  set lrcmdiv [expr {1/(1 + [lrcmd $lrcdVec1 $lrcdVec2 $infoV $loSt])}]
  if {$out} {
    puts $loSt "\nlig-rec contact distance Manhattan distance index: $lrcmdiv"
    }
  return $lrcmdiv
  }   ;# lrcmdi


#|-proc lrcp {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} :
#|  -calculates the LRCP parameter from two L-R contact vectors .
#|  -definitions :
#|    -LRCP :
#|      -ligand-receptor contact projection or the dot product of two
#|       _ normalized lrcd vectors ;;
#|  -arguments :
#|    -lrcdVec1 :-L-R contact vector for reference complex .
#|      -list formatted to build an array as returned by lrcdVec_sum ;
#|    -lrcdVec2 :-L-R contact vector for the complex to compare to ref .
#|      -list formatted to build an array as returned by lrcdVec_sum ;
#|    -infoV :-information string about the systems .
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;
#|        -'none':-will avoid any output ;;;;;
proc lrcp {lrcdVec1 lrcdVec2 {infoV ""} {loSt stdout}} {
  set lrcdVec1 [lrcdVec_norm $lrcdVec1 "none"]
  set lrcdVec2 [lrcdVec_norm $lrcdVec2 "none"]
  array unset lrcdArr1
  array unset lrcdArr2
  array set lrcdArr1 $lrcdVec1
  array set lrcdArr2 $lrcdVec2
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {[llength [array names lrcdArr1]] != [llength [array names lrcdArr2]]} {
    if {$out} {puts $loSt "lrcp: Error: Incompatible LRCD vectors"}
    return "-"
    }
  if {$out} {
    puts $loSt "\n+++ligand-receptor contact projection (LRCP) of two complexes +++"
    puts $loSt "\nuser info: $infoV"
    puts $loSt "types considered in the L-R complex 1: [array names lrcdArr1]"
    puts $loSt "types considered in the L-R complex 2: [array names lrcdArr2]"
    puts $loSt "\n component \t complex 1 \t complex 2 \t cmplx 2 · cmplx 1"
    }
  set sum 0.0
  foreach type [array names lrcdArr1] {
    set prod [expr {$lrcdArr2($type)*$lrcdArr1($type)}]
    if {$out} {puts $loSt "$type\t$lrcdArr1($type)\t$lrcdArr2($type)\t$prod"}
    set sum [expr {$sum + $prod}]
    }
  set lrcpv [expr {$sum}]
  if {$out} {
    puts $loSt "\nlig-rec contacts projection: $lrcpv"
    puts $loSt "\nlrcp: Done."
    }
  return $lrcpv
  }   ;# lrcp


#|-proc lrcd_pwMat {selIdL cutoff {outPref "lrcd_pwMat"} {vecDist nd}
#|                              _ {lrcdMatL {dct cmd lapa}} {loSt stdout}} :
#|  -generates a matrix with pairwise values of LRCD for a set of complexes .
#|  -allows selecting different variants of the LRCD parameter
#|   _ (through the vecDist argument) .
#|  -the output of the procedures used to calculate the matrix is log out into
#|   _ a text file .
#|  -arguments :
#|    -selIdL :
#|      -list of selId's for the selInfo array .
#|      -each selId will require the title, ligId, ligSelTxt, recId and
#|       _ recSelTxt keys declared in the selInfo array :
#|        -example: selInfo(selId,ligSelTxt) ;;
#|    -cutoff :
#|      -maximum distance limit (in Angstroms) to consider contacts and
#|       _ to scale some parameters .
#|      -recommended value :-4.5 ;;
#|    -outPref :
#|      -prefix for the output log file name .
#|      -the output log file will contain the output of several procedures
#|       _ called during the lrcd calculation .
#|      -default value :-'lrcd_pwMat' ;
#|      -acceptable values :
#|        -'stdout' :-output to the console, no log file is created ;;;
#|    -vecDist :-lrc vector distance type .
#|      -choose the method to calculate the difference of the two lrc vectors .
#|      -acceptable values :
#|        -"d" :-Cartesian distance between the lrc vectors ;
#|        -"nd" :-Cartesian distance between the normalized vectors ;
#|        -"di" :-Cartesian distance index [= 1/(1 + LRCND)] ;
#|        -"md" :-Manhattan distance between the normalized vectors ;
#|        -"mdi" :-Manhattan distance index [= 1/(1 + LRCMD)] ;
#|        -"p" :-dot product (internal projection) of the normalized vectors ;
#|        -other methods may be implemented in future versions ;
#|      -default value :-"nd" ;;
#|    -lrcdMatL :
#|      -list of flags indicating which lrcdMat procedures are to be
#|       _ considered in the creation of the lrcd matrices .
#|      -accepatble values :
#|        -a tcl list with one or more of the flags 'cdt', 'cmd' and 'lapa' .
#|        -other procedures may be incorporated in future versions ;
#|      -default value :-{dct cmd lapa} ;;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;;;
#|  -notes :
#|    -requires the global selInfo array to be propertly set .
#|    -in this case only the lrcdVec of type 'sum' is considered, but in
#|     _ future version should be able to choose other types of vectors .
#|    -other types of components (keys) in the lrcdMat should be incorporated
#|     _ in future versions .
#|    -currently a square symmetric matrix is calculated, but it may be changed
#|     _ to calculate on the upper or lower diagonal matrix ;;
proc lrcd_pwMat {selIdL cutoff {vecDist nd} {outPref "lrcnd_pwMat"} \
                               {lrcdMatL {dct cmd lapa}} {loSt stdout}} {
  global selInfo lrcd_version
  array unset lrcdVecs
  if {$outPref == "lrcnd_pwMat"} {set outPref "lrc${vecDist}_pwMat"}
  if {($vecDist != "d") && ($vecDist != "nd") && ($vecDist != "di") \
      && ($vecDist != "md") && ($vecDist != "mdi") && ($vecDist != "p")} {
    puts "lrcd_pwMat: Error! Acceptable values of vecDist are d, nd, di, md, mdi, and p."
    return
    }
  if {$outPref == "stdout"} {set log stdout} else {
    set outPref "${outPref}_$cutoff-[join $lrcdMatL "-"].log"
    set log [open $outPref w]
    }
  puts $loSt "\n+++ pairwise LRCD Matrix for a set of complexes +++"
  puts $loSt "\nselIds (complexes): $selIdL"
  puts $loSt "flags for LRCD matrices: $lrcdMatL"
  puts $loSt "cutoff: $cutoff"
  puts $loSt "output log: $outPref\n"
# calculates the lrcdVec for each complex (selId list)
  puts $log "\n+++ pairwise LRCD Matrix for a set of complexes +++"
  puts $log "lrcd.tcl v-$lrcd_version, lrcd_pwMat, [exec date]"
  foreach selId $selIdL {
    set lrcdMat {}
    puts $log "\nCalculating LRCD matrix for complex specified in $selId"
    foreach flag $lrcdMatL {
      set lrcdMat [lrcdMat_$flag $lrcdMat $selId $cutoff $log]
      }
    puts $log "\nCalculating LRCD vector for complex specified in $selId"
    set lrcdVecs($selId) [lrcdVec_sum $lrcdMat $selInfo($selId,title) $log]
    }
# calculates the lrcd matrix for each pair of complexes
  foreach selId $selIdL {
    puts -nonewline $loSt "\t$selId"
    }
  puts $loSt ""
  foreach selIdi $selIdL {
    puts -nonewline $loSt "$selIdi"
    foreach selIdj $selIdL {
      puts $log "\nCalculating LRCD between complexes in $selIdi and $selIdj"
      set lrcd_ij [lrc$vecDist $lrcdVecs($selIdi) $lrcdVecs($selIdj) \
                   "$selInfo($selIdi,title) vs $selInfo($selIdj,title)" $log]
      if {$lrcd_ij == "-"} {
        puts -nonewline $loSt "\t---"
      } else {
        puts -nonewline $loSt "[format "\t%.3f" $lrcd_ij]"
        }
      }
    puts $loSt ""
    }
  if {$outPref != "stdout"} {close $log}
  puts "\nlrcd_pwMat: Done."
  }   ;# *** lrcd_pwMat ***

#|-proc lrcd_dlgTab {selIdL selIdRef cutoff {vecDist nd}
#|       _ {outPref "lrcnd_dlgTab"} {lrcdMatL {dct cmd lapa}} {loSt stdout}} :
#|  -generates a table with lrcd values for a set of complexes contained in a
#|   _ AutoDock (AD) .dlg file against a reference complex .
#|  -note that AutoDock4 (AD4) .dlg files may also be processed provided that
#|   _ no flexible residues are considered in the AutoDock calculation .
#|  -each row of the table corresponds to different ligand structures derived
#|   _ from the cluster analysis performed by AD (best structure of each
#|   _ cluster) .
#|  -allows selecting different variants of the LRCD parameter
#|   _ (through the vecDist argument) .
#|  -arguments :
#|    -selIdL :
#|      -list of selIds for the selInfo array refering to lig-rec complexes .
#|      -each selId is supposed to contain several frames with different
#|       _ conformations for the complex ;
#|    -selIdRef :
#|      -selId for the selInfo array refering to a lig-rec complex which is
#|       _ considered as a common reference for the complexes in selIdL ;
#|    -cutoff :
#|      -maximum distance limit (in Angstroms) to consider contacts and
#|       _ to scale some parameters .
#|      -recommended value :-4.5 ;;
#|    -outPref :
#|      -prefix for the output log file name .
#|      -the output log file will contain the output of several procedures
#|       _ called during the lrcd calculation .
#|      -default value :-'lrcnd_dlgTab' ;
#|      -acceptable values :
#|        -'stdout' :-output to the console, no log file is created ;;;
#|    -vecDist :-lrc vector distance type .
#|      -choose the method to calculate the difference of the two lrc vectors .
#|      -acceptable values :
#|        -"d" :-Cartesian distance between the vectors ;
#|        -"nd" :-Cartesian distance between the normalized vectors ;
#|        -"di" :-Cartesian distance index [= 1/(1 + LRCND)] ;
#|        -"md" :-Manhattan distance between the normalized vectors ;
#|        -"mdi" :-Manhattan distance index [= 1/(1 + LRCMD)] ;
#|        -"p" :-dot product (internal projection) of the normalized vectors ;
#|        -other methods may be implemented in future versions ;;
#|    -lrcdMatL :
#|      -list of flags indicating which lrcdMat procedures are to be
#|       _ considered in the creation of the lrcd matrices .
#|      -accepatble values :
#|        -a tcl list with one or more of the flags 'cdt', 'cmd' and 'lapa' ;
#|      -default value :-{dct cmd lapa} ;;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;;;
#|  -notes :
#|    -the matrix reported contains in the first column the selId, so each row
#|     _ lists the LRCD for all the conformations for each comples, thus the
#|     _ table should be transposed to have each complex on each column ;;
proc lrcd_dlgTab {selIdL selIdRef cutoff {vecDist nd} {outPref "lrcnd_dlgTab"} \
                  {lrcdMatL {dct cmd lapa}} {loSt stdout}} {
  global selInfo lrcd_version
  if {$outPref == "lrcnd_dlgTab"} {set outPref "lrc${vecDist}_dlgTab"}
  if {($vecDist != "d") && ($vecDist != "nd") && ($vecDist != "di") \
      && ($vecDist != "md") && ($vecDist != "mdi") && ($vecDist != "p")} {
    puts "lrcd_dlgTab: Error! Acceptable values of vecDist are d, nd, di, md, mdi, and p."
    return
    }
  if {$outPref == "stdout"} {set log stdout} else {
    set outPref "${outPref}_$selIdRef-$cutoff-[join $lrcdMatL "-"].log"
    set log [open $outPref w]
    }
  puts $loSt "\n+++ LRCD Table for a set of Autodock .dlg docking results +++"
  puts $loSt "\nselIds (complexes): $selIdL"
  puts $loSt "ref complex: selId: $selIdRef; title: $selInfo($selIdRef,title)"
  puts $loSt "flags for LRCD matrices for each complex: $lrcdMatL"
  puts $loSt "cutoff: $cutoff"
  puts $loSt "output log: $outPref\n"
# calculate the lrcdVec for the reference complex
  set lrcdMat {}
  puts $log "\n+++ LRCD Table for a set of Autodock .dlg docking results +++"
  puts $log "lrcd.tcl v-$lrcd_version, lrcd_dlgTab, [exec date]"
  puts $log "\nCalculating LRCD matrix and vector for ref complex: $selIdRef"
  foreach flag $lrcdMatL {
    set lrcdMat [lrcdMat_$flag $lrcdMat $selIdRef $cutoff $log]}
  set lrcdVecRef [lrcdVec_sum $lrcdMat "ref: $selInfo($selIdRef,title)" $log]
# calculate the lrcdVec for each complex and frame and the lrcd values
  foreach selId $selIdL {
    set id $selInfo($selId,ligId)
    puts $log "\nCalculating LRCD for each structure in complex $selId"
    puts -nonewline $loSt "$id\t$selId\t"
    for {set f 0} {$f < [molinfo $id get numframes]} {incr f} {
      mol top $id
      animate goto $f
      set lrcdMat {}
      puts $log "\nCalculating $selId frame: $f"
      foreach flag $lrcdMatL {
        set lrcdMat [lrcdMat_$flag $lrcdMat $selId $cutoff $log]}
      set lrcdVec [lrcdVec_sum $lrcdMat "$selInfo($selId,title)" $log]
      set lrcdi [lrc$vecDist $lrcdVec $lrcdVecRef "$selInfo($selId,title) (frm $f) vs $selInfo($selIdRef,title)" $log]
      if {$lrcdi == "-"} {
        puts -nonewline "---\t"
      } else {
        puts -nonewline [format "%.2f\t" $lrcdi]
# storing the calculated lrcd in the user field of each frame
        [atomselect $id all frame $f] set user [format "%.2f\t" $lrcdi]
        }
      }
    animate goto 0
    puts $loSt ""
    }
# this updates the ranges of the user field from 0.0 to $max
  foreach selId $selIdL {
    mol scaleminmax $id 1 
    }
  if {$outPref != "stdout"} {close $log}
  puts "\nlrcd_dlgTab: Done."
  }   ;# *** -lrcd_dlgTab ***



#|-proc lrcd_dbTab {selIdL lrcdVecArr {cutoff 4.5} {vecDist d}
#|                _ {outPref "lrcdi_dbTab"} {lrcdMatL {dct cmd lapa}}
#|                _ {loSt stdout}} :
#|  -calculates the LRCD parameter or one of its variants for a complex or
#|   _ a set of docking poses against a LRCD vector data base of PDB complexes .
#|  -arguments :
#|    -selIdL :
#|      -list of selIds for the selInfo array refering to lig-rec complexes .
#|      -each selId may contain several frames with different
#|       _ conformations for the complex ;
#|    -lcdVecArr :
#|      -Tcl array with a set of lrcdVec arrays pre-calculated for a set of
#|       _ complexes .
#|      -it may be created declaring an array and using the
#|       _ lrcdVec_sum procedure or may be read from a file using the
#|       _ lrcdVec_readDB procedure ;
#|    -cutoff :
#|      -maximum distance limit (in Angstroms) to consider contacts and
#|       _ to scale some lrcd parameters .
#|      -default value :-4.5 ;;
#|    -vecDist :-lrc vector distance type .
#|      -choose the method to calculate the difference of the two lrc vectors .
#|      -acceptable values :
#|        -"d" :-Cartesian distance between the vectors (LRCD) ;
#|        -"nd" :-Cartesian distance between the normalized vectors (LRCND) ;
#|        -"di" :-Cartesian distance index [LRCDI = 1/(1 + LRCND)] ;
#|        -"md" :-Manhattan distance between the normalized vectors (LRCMD) ;
#|        -"mdi" :-Manhattan distance index [LRCMDI= 1/(1 + LRCMD)] ;
#|        -"p" :-dot product (projection) of the normalized vectors (LRCP) ;
#|        -other methods may be implemented in future versions ;;
#|    -outPref :
#|      -prefix for the output log file name .
#|      -the output log file will contain the output of several procedures
#|       _ called during the lrcd calculation .
#|      -default value :-'lrcdi_dbTab' ;
#|      -acceptable values :
#|        -'stdout' :-output to the console, no log file is created ;;;
#|    -lrcdMatL :
#|      -list of flags indicating which lrcdMat procedures are to be
#|       _ considered in the creation of the lrcd matrices .
#|      -accepatble values :
#|        -a tcl list with one or more of the flags 'cdt', 'cmd' and 'lapa' ;
#|      -default value :-{dct cmd lapa} ;;
#|    -loSt :-output stream (channel) for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id refering to an already open output text file ;;;;
proc lrcd_dbTab {selIdL lrcdVecArr {cutoff 4.5} {vecDist d} \
                 {outPref "lrcdi_dbTab"} {lrcdMatL {dct cmd lapa}} \
                 {loSt stdout}} {
# initializing and input check
  global selInfo lrcd_version
  array unset lrcdVecL
  if {$outPref == "lrcndi_dbTab"} {set outPref "lrc${vecDist}_dbTab"}
  if {($vecDist != "d") && ($vecDist != "nd") && ($vecDist != "di") \
      && ($vecDist != "md") && ($vecDist != "mdi") && ($vecDist != "p")} {
    puts "lrcd_dbTab: Error! Acceptable values of vecDist are d, nd, di, md, mdi, and p."
    return
    }
  if {$outPref == "stdout"} {set log stdout} else {
    set outPref "${outPref}_[join $selIdL "-"]-$cutoff.log"
    set log [open $outPref w]
    }
  array set lrcdVecDB $lrcdVecArr
# printing info
  puts $loSt "\n+++ LRCD Table for a set of complexes against a lrcdVec data base +++"
  puts $loSt "\nselIds (complexes): $selIdL"
  puts $loSt "\nIDs for the complexes DB: [array names lrcdVecDB]"
  puts $loSt "flags for LRCD matrices for each complex: $lrcdMatL"
  puts $loSt "LRCD parameter type: lrc$vecDist"
  puts $loSt "cutoff: $cutoff"
  puts $loSt "output log: $outPref"
  puts $log "\n+++ LRCD Table for a set of complexes against a lrcdVec data base +++"
  puts $log "lrcd.tcl v-$lrcd_version, lrcd_dbTab, [exec date]"
  puts $log "\nCalculating LRCD table for the complexes: $selIdL"
# creating header for output table in the console
  puts -nonewline $loSt "\nNo.\tRef-PDBID-LIG\t "
  foreach selId $selIdL {
    puts -nonewline "$selId"
    set ligId $selInfo($selId,ligId)
    for {set f 0} {$f < [molinfo $ligId get numframes]} {incr f} {
      puts -nonewline ":$f\t"
      }
    }
  puts $loSt ""
# calculating lrcdVec for the selIds in selIdL
  foreach selId $selIdL {
    set ligId $selInfo($selId,ligId)
    puts $log "Calculating lrcdVec for complexes in selId: $selId"
    for {set f 0} {$f < [molinfo $ligId get numframes]} {incr f} {
      mol top $ligId
      animate goto $f
      set lrcdMat {}
      foreach flag $lrcdMatL {
        set lrcdMat [lrcdMat_$flag $lrcdMat $selId $cutoff $log]}
      set lrcdVecL(${selId}_$f) [lrcdVec_sum $lrcdMat "$selInfo($selId,title)" $log]
      }
    animate goto 0
    }
# caclulating lrcd for each frame in each $selId taking each ref vector
  set ic 0
  foreach refCompId [array names lrcdVecDB] {   ;# running over ref complexes
    puts $log "\nEvaluating LRCD with reference lrcdVec: $refCompId"
    set lrcdVecRef $lrcdVecDB($refCompId)
    puts -nonewline $loSt "$ic\t$refCompId\t"
    foreach selId $selIdL {
      set ligId $selInfo($selId,ligId)
      puts $log "Calculating lrcd for complexes in selId: $selId"
      for {set f 0} {$f < [molinfo $ligId get numframes]} {incr f} {
        puts $log "\nCalculating $selId frame: $f"
        set lrcdi [lrc$vecDist $lrcdVecL(${selId}_$f) $lrcdVecRef "$selInfo($selId,title) (frm $f) vs ref complex $refCompId" $log]
        if {$lrcdi == "-"} {
          puts -nonewline "---\t"
        } else {
          puts -nonewline [format "%.3f\t" $lrcdi]
          }
        }
      }
    puts $loSt ""
    incr ic
    }
  if {$outPref != "stdout"} {close $log}
  puts "\nlrcd_dbTab: Done."
  }   ;# lrcd_dbTab



#|-proc lrcdVec_len {lrVec {loSt stdout}} :
#|  -returns the norm of a vector of measurements of ligand-receptor binding .
#|  -arguments :
#|    -lrVec :-list formated to build an array with L-R measurements as
#|     _ returned by ligRecContVec_SD ;
#|    -loSt :-output stream for log messages ;;;
proc lrcdVec_len {lrVec {loSt stdout}} {
  array unset lrcv
  array set lrcv $lrVec
  set types [array names lrcv]
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {$out} {
    puts $loSt "\n+++ Norm of a vector describing ligand-receptor contacts +++"
    puts $loSt "types (components) in the vector: $types"
    }
  set sum 0.0
  foreach t $types {
    set sum [expr {$sum + $lrcv($t)*$lrcv($t)}]
    }
  set norm [expr {sqrt($sum)}]
  if {$out} {
    puts $loSt "\nvector norm: $norm"
    puts $loSt "\nlrcdVec_len: Done."
    }
  return $norm
  }   ;# lrcdVec_len

#|-proc lrcdVec_norm {lrcdVec {loSt stdout}} :
#|  -returns a lrcd vector normalized .
#|  -arguments :
#|    -lrcdVec :
#|      -list formatted to build an array with preprocessed L-R contacts
#|       _ measurements (a type of contact per component) ;
#|    -loSt :-output stream for log messages .
#|      -default value :-'stdout' ;
#|      -acceptable values :
#|        -a channel Id referring to an already open output text file .
#|        -'none':-will avoid any output ;;;;;
proc lrcdVec_norm {lrcdVec {loSt stdout}} {
  array unset lrcdv
  array set lrcdv $lrcdVec
  set types [array names lrcdv]
  if {$loSt == "none"} {set out 0} else {set out 1}
  if {$out} {
    puts $loSt "\n+++ ligand-receptor contact distance vector normalization +++"
    puts $loSt "categories in the lrcd vector: $types"
    puts $loSt "\n type \t dist sum norm"
    }
# calculate the magnitude of the vector
  set norm [lrcdVec_len $lrcdVec $loSt]
  foreach type $types {
    set lrcdv($type) [expr {$lrcdv($type)/$norm}]
    if {$out} {puts $loSt "$type\t$lrcdv($type)"}
    }
  if {$out} {puts $loSt "\nlrcdVec_norm: Done."}
  return [array get lrcdv]
  }   ;# lrcdVec_norm


