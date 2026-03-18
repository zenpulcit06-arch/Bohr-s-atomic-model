#this is shell script run.sh
#!/bin/bash
gcc Simulation.c -o simulation
echo "Enter orbit of electron"
read orbit
echo "$orbit" | ./simulation
python Animation.py
