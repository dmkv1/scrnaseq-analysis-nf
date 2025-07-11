#!/bin/bash
set -e

nohup ./inferCNV.R > inferCNV.log 2>&1 & disown
echo $! > inferCNV.pid
