#/bin/bash

nohup ./odt_multi_realz.sh "steady_state" "nml_mphy_inj0.3_aero"  "inj0.3_aero"   &> bat_inj0.3.log &
nohup ./odt_multi_realz.sh "steady_state" "nml_mphy_inj0.5_aero"  "inj0.5_aero"   &> bat_inj0.5.log &
nohup ./odt_multi_realz.sh "steady_state" "nml_mphy_inj0.8_aero"  "inj0.8_aero"   &> bat_inj0.8.log &
nohup ./odt_multi_realz.sh "steady_state" "nml_mphy_inj1_aero"    "inj1.0_aero"   &> bat_inj1.log &
nohup ./odt_multi_realz.sh "steady_state" "nml_mphy_inj2_aero"    "inj2.0_aero"   &> bat_inj2.log &
nohup ./odt_multi_realz.sh "steady_state" "nml_mphy_inj3_aero"    "inj3.0_aero"   &> bat_inj3.log &
nohup ./odt_multi_realz.sh "steady_state" "nml_mphy_inj6_aero"    "inj6.0_aero"   &> bat_inj6.log &
nohup ./odt_multi_realz.sh "steady_state" "nml_mphy_inj10_aero"   "inj10.0_aero"  &> bat_inj10.log &

