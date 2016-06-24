all:
	f2py -m opac -c opacity_routine.f wrapper.f95

get:
	-wget http://www2.mpia-hd.mpg.de/home/henning/Dust_opacities/Opacities/Code/opacity.zip && unzip -o opacity.zip && rm opacity.f~ opacity.zip 
	
clean:
	-rm semenov_opacs.pyc
	-rm kR.out
	-rm opacity.f
	-rm opacity.inp
	-rm opacity_routine.f
	-rm Readme.txt

clobber: clean
	-rm -rf opac.so opac.so.dSYM
	-rm kP_h2001.dat
	-rm kR_h2001.dat
