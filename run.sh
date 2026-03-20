#this is shell script run.sh
#!/bin/bash
rm -rf simulation
gcc -I. sim_essen/acc_cal/acc.c sim_essen/verlet/integrator.c Simulation.c -o  simulation -lm
if [ $? -ne 0 ]; then
  echo "Error" >&2
  exit 1
fi
echo "Enter orbit of electron"
read orbit
echo "$orbit" | ./simulation
python Animation.py
