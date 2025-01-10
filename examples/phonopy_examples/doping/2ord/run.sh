for i in {001..062}
do
ln -s /home/jxchen/Data/vasp/dopgraphene/perfect/CHGCAR $i/CHGCAR
ln -s /home/jxchen/Data/vasp/dopgraphene/perfect/WAVECAR $i/WAVECAR
	cd $i
	mpirun -np 48 vasp_std > job.log 
	cd ..
done
