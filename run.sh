#this is shell script run.sh
#!/bin/bash
rm -f simulation
rm orbit.csv
rm photon.csv

gcc -I. sim_essen/acc_cal/acc.c sim_essen/verlet/integrator.c Simulation.c -o  simulation -lm

if [ $? -ne 0 ]; then
  echo "Error" >&2
  exit 1
fi

echo "Enter orbit of electron"
read orbit
echo "$orbit" | ./simulation

if [[ ! -f "orbit.csv" ]]; then
  echo "Simulation failed"
  exit 1
fi
python Animation.py

rm photon.csv
rm orbit.csv
rm -rf simulation
