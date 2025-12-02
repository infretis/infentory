#!/bin/bash
clear
rm -rf infretis_data*.txt amsworker* *load/ worker* run?/ infretis_?.toml pattern.txt sim.log infretis_init.log
cp infretis.toml_bkp fresh_start.toml
inft infinit -toml fresh_start.toml
