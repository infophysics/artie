#
# Source this setup from within your working directory
# where you plan to compile a stand-alone GEANT application
#
export G4WORKDIR=$PWD
source /dune/app/users/mulhearn/standalone/geant4.10.05.p01-install/geant4make.sh 
source /dune/app/users/mulhearn/standalone/root-6.18.00-build/bin/thisroot.sh

#
# for neutronHP
#
unset G4NEUTRONHP_SKIP_MISSING_ISOTOPES 
unset G4NEUTRONHP_DO_NOT_ADJUST_FINAL_STATE
unset G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION
unset G4NEUTRONHP_NELECT_DOPPLER
unset G4NEUTRONHP_PRODUCE_FISSION_FRAGMENTS

#source G4NEUTRONHPDATA=/home/yashwanth/G4NDL4.5

#### G4NEUTRONHP_SKIP_MISSING_ISOTOPES=1
#### export G4NEUTRONHP_SKIP_MISSING_ISOTOPES

#### G4NEUTRONHP_DO_NOT_ADJUST_FINAL_STATE=1
#### export G4NEUTRONHP_DO_NOT_ADJUST_FINAL_STATE
   
#### G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1
#### export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION
  
#### G4NEUTRONHP_NELECT_DOPPLER=1
#### export G4NEUTRONHP_NELECT_DOPPLER
  
#### G4NEUTRONHP_PRODUCE_FISSION_FRAGMENTS=1
#### export G4NEUTRONHP_PRODUCE_FISSION_FRAGMENTS
#
# for Bertini cascade
#
unset G4CASCADE_USE_PRECOMPOUND  
#### G4CASCADE_USE_PRECOMPOUND=1
#### export G4CASCADE_USE_PRECOMPOUND             

env |grep G4

