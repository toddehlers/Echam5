#! /bin/csh -fx
#PBS -l nodes=2:ppn=8       
#PBS -l walltime=900:00:00 
##PBS -l mem=1000mb     
#PBS -j oe       
#PBS -N S2.3
#PBS -M rfeng@umich.edu 
################################################################################
# Job file to run echam model on toucan 
# ======================================
setenv EXP S2.3           # experiment identifier
set RES = 106             # resolution 
set LEV = 19		 # how many levels you like
set RERUN = .true.       # for initial run, set to .false. for restart run set to .true.
set REMON = 6            # only has effect on restart run, specify the month you'd like to restart the model. This variable should match the month of restart file +1
set YRST = 1991          # the year to start/restart the run
set YRED = 2020          # the year to stop the run
setenv YEAR 1998
set lmlo = .false.       # whether to run with slab ocean, currently unavailable
set lamip = .false.      # whether to run with cyclic SST, AMIP mode, transient climate run 
#
set UT  = /condor/data18/rfeng/echam5-wiso.v1.1         #home directory of executable
set MODEL = /condor/data18/rfeng/echam5-wiso.v1.1/bin/echam5  			   #exact location of the executable
set UTF = /condor/data18/rfeng/echam5-wiso.v1.1/run     #upper level run directory where all the runs are stored
set DAT = $UTF/Eocene/${EXP}				   #run directory where all the output are found
setenv DPATH ${UTF}/Eocene/${EXP}/			   #path of run directory
set DRERUN = $DAT/rerun 		           #where restart file for each year are stored, in case of branching a run 
set MPIROOT = /opt/packages/openmpi/1.4.4/intel/mx  #path of MPI root
set MPI_BIN = $MPIROOT/bin			   #MPI library
#################################################################################directory of boundary dataset
#
set INI = /condor/data18/rfeng/Data_echam/T${RES}			   #initial data for atmosphere
set INIAMIP = ${INI}/amip2          #initial data for sst
set INIWISO = /condor/data18/rfeng/Data_echam/WISO/T${RES}           #initial data for water isotope
#end of user input
#################################################################################
#
set NCPUS = 16 					   #how many cpus are used in this run, equals nodesxppn
set NPROCA = 4					   #cpus to process task A
set NPROCB = 4					   #cpus to process task B, which is defined by domain decompostion during parallel computing.  
set NTHREADS = 1				   #Only this option available, it seems unable to do a hybrid run
set NPROMA = 64 				   #vector length during computing
#
#################################################################################
#
 setenv ECHAM5_THREADS $NTHREADS
 setenv  OMP_NUM_THREADS $NTHREADS 
 unlimit
 setenv KMP_MONITOR_STACKSIZE 1m
 setenv KMP_STACKSIZE 32m
 setenv OMP_DYNAMIC FALSE
#
################################################################################
#
if ( ! -e $DPATH ) then
   mkdir $DPATH
endif

if ( ! -e $DAT ) then
   mkdir $DAT
endif

if ( ! -e $DRERUN ) then
   mkdir $DRERUN
endif
#################################################################################
#
cd $DPATH                          # output and rerun files are written into $DPATH
#
################################################################################
#for initialization of the atmosphere and land surface
rm -f unit.?? sst* ice* rrtadata *.codes atmout
#
#############################################the following file name doesn't matter, only unit.XX matters, make sure not change the ".XX" if you want to read in certain variable sets
#
ln -s  ${INI}/T${RES}L${LEV}_jan_spec.nc          unit.23
ln -s  ./T${RES}_S2.3.nc                 unit.24  # this directory should be changed if any change is made to use other surface boundary dataset
#
################################################################################
#for initialization of land surface
ln -s  ${INI}/T${RES}_O3clim2.nc    unit.21
ln -s  ./T${RES}_VLTCLIM.nc    unit.90
ln -s  ./T${RES}_VGRATCLIM.nc  unit.91
ln -s  ./T${RES}_TSLCLIM2.nc   unit.92
ln -s  ${INI}/surrta_data       rrtadata
#
################################################################################
## for climatological sst and ice (LAMIP=F) use:
if ( $lamip == .false. )then
ln -s  ./T${RES}_amip2sst_clim.nc    unit.20
ln -s  ./T${RES}_amip2sic_clim.nc    unit.96
endif
## for slab ocean module, currently unavailable
if ( $lmlo == .true. )then
ln -s  ./T${RES}_heatflux_test.nc 	      unit.42
endif
##
################################################################################
# for AMIP (variable) sst and ice (LAMIP=T) use:
if ( $lamip == .true. )then
set i = $YRST
@ i = $i - 1
while ($i <= $YRED)
ln -s  ${INIAMIP}/T${RES}_amip2sst_${i}.nc  sst${i}
ln -s  ${INIAMIP}/T${RES}_amip2sic_${i}.nc  ice${i}
@ i++
end
endif
#
################################################################################
#need for initialization of water isotope
ln -s  ./T${RES}_wisosw_d.nc        unit.25
#
################################################################################
#  namelist control variables and output control for grid variables
#  spectral variables are written out by default except liq. water
#  for production runs set LABORT=.FALSE.
#
while ( $YEAR <= $YRED )
#######################################################
set namelist = namelist.echam
@ ied = $YEAR + 1
echo "cat >! ${namelist} << EOF"
####################################################cycling for a year
###########################only need for the first time to start from months other than Jan
if ( $RERUN == .true. ) then
     set mon = $REMON 
endif
###########################
# only run for one month
while ( $mon < 13 )
if ( $mon < 10 ) then
     set MON = "0${mon}"
else
     set MON = "${mon}"
endif
#####################################################
#namelist of ECHAM5-wiso
cat >! ${namelist} << EOF
&RUNCTL
  LRESUME=$RERUN,
 out_datapath = "$DPATH"
 out_expname  = "$EXP"
  DT_START  = $YRST,01,01,0,0,0                 
  DT_STOP   = $YRED,01,01,0,0,0
  out_filetype = 2  
  PUTDATA   = 6,'hours','last',0                             
  TRIGJOB(1)= 1,'months','last',0                          
  NSUB=0,						  
  putrerun  = 1, 'months', 'last', 0			 
  LAMIP=$lamip,						
  LMLO=$lmlo,					
  LABORT=F,					
  NPROCA=${NPROCA}
  NPROCB=${NPROCB}
  NPROMA=${NPROMA}
  DELTA_TIME=360.,
/
&RADCTL
  CO2VMR=1120.E-06,
/
&WISOCTL
  NWISO=3,
/
EOF
#
###################################################################
#inherent ECHAM5 parallel computing environment variable, need to be the same as environmental nthreads
setenv ECHAM5_THREADS $NTHREADS
#mpirun -np $NCPUS --mca btl mx,sm,self $MODEL >> fort.6_${YEAR}${MON} || echo "run failed" && exit 1
mpirun -np $NCPUS $MODEL >> fort.6_${YEAR}${MON} || echo "run failed" && exit 1
#$MPI_BIN/mpiexec --prefix $MPIROOT -np $NCPUS $MODEL >>& fort.6_${YEAR}${MON}
###########################################end of the model run for one year
#
set statenv = $status

if ( $statenv != 0) then
    exit 
endif

echo "end of run, status is ${statenv}" >>& fort.6_${YEAR}${MON}
cp rerun_${EXP}_echam ${DRERUN}/rerun_${EXP}_echam_${YEAR}${MON}
set RERUN = .true.
#switch the rerun option to do continuous run
##########################################
#
@ mon = $mon + 1
end
#
##########################################start postprocessing
# variables used in subjob2
#set YEARED = $YEAR
#set DAY = 01
#cp ${UTF}/subjob2 .    			#copy subjob2 to run directory and submit subjob2 to do postprocessing
#qsub -V ./subjob2
##########################################
#
@ YEAR = $YEAR + 1		        #start the run for another year
set mon  = 1
#########################################end of the whole run
end
################################################################################
