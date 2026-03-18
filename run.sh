#this is shell script run.sh
gcc Simulation.c -o simulation 
echo "1" | ./simulation
python Animation.py
