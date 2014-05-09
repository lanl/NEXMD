##SCRIPT TO STEP THROUGH VALUES OF Ex in NAESMD INPUT AND SORT EXCITATION ENERGIES INTO OUTPUT FILES
#!/bin/bash
rm *.out
li=$2
nx=$3
stepx=$1
#echo 'Step size for R (0):'
#read stepx
#echo 'Start/Final steps (0 1):'
#read li nx
#echo 'Got it'
nstates=$(grep -oP 'n_exc_states_propagate=\K[0-9]*' input.ceon) ##Get the number of states being calculated
echo "Doing $nstates states"
nstates=$((nstates+1))
for ((l=li;l<$((nx+1));l++))
do
##Replace value in input file and run NAESMD
x=$(echo "scale=10; $l*$stepx " | bc -l)
R=$(echo "scale=10;  e(l(0.5*$x+1)-l(1-$x))" | bc -l) ##Calculate the new Ex values


echo "step $x E_$x.out Epsilon=$R"

sed "s\ceps=[\.0-9]*\ceps=$R\g" -i input.ceon ##Replace Ex value in input file 
./sqmceonaesmd.exe > "E_$x.out" ##Run NAESMD

##Comb output for ground and excitation energies and print to files
Egr=$(grep -A 1 'Total energy of the ground state' E_$x.out | grep -oP '0[ ]+[-0-9\.]+'| awk '{print $2}') #grep gs total energy
echo "Egr=$Egr"
echo "$x $Egr" >> Egr.out #write to file
for ((n=1;n < $((nstates)) ;n++))
do
E=$(grep -A $nstates 'Transition Dipole' E_$x.out | grep -oP -m 1 " $n[ ]+\K\-*[0-9]\.*[0-9]*") #grep the excitation energies out
echo "$x $E" >> Eex_$n.out #write them to files
echo "E$n=$E"
done

echo "Current R= $x"
done
echo "Finished"
