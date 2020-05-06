#!/usr/bin/env bash

odir=./test
winds="0.5 3 8 15" # 7.2 14.2"
wazis=$(seq 0 30 359)
szas=3 #$(seq 0 20 60)
slope="" #--slope

for wind in $winds; do
 for wazi in $wazis;do
  for sza in $szas;do

    echo $sza $wind $wazi
    odir_=$odir/sunglint/$sza
    mkdir -p $odir_

    filename=$odir_/coxmunk_fig_sza"$sza"_ws"$wind"_wazi"$wazi"_bh2006$slope.png

    coxmunk $sza $wind --stats bh2006 --wind_azi $wazi --shadow $slope --figname $filename

  done
 done
done