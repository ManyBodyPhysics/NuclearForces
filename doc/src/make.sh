cwd=$(pwd)
cd obemodels && bash make.sh obemodels
cd ${cwd}
cd forces && bash make.sh forces
cd ${cwd}
cd eft && bash make.sh eft
cd ${cwd}
cd phenomenologyforces && bash make.sh phenomenologyforces
cd ${cwd}
cd scatteringtheory && bash make.sh scatteringtheory
cd ${cwd}
cd introduction && bash make.sh introduction
cd ${cwd}
