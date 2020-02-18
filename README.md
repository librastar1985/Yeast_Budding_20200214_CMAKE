# budding_cmake
After forking the code, put on crc. then: navigate to the folder with build scripts. run ./build_script.sh <your_build_folder_here> I tend to use ./build_script.sh build The executable and static libraries are in "build" and code is in source.

When pushing back to github, remove your build folder using rm -rf <your_build_folder_here> since you only should push the source code and make files. Always push while in "budding_cmake"

Place the budding_cmake folder in a clean folder under your CRC directory.
Create a folder called Simulation.
    -> Under Simulation, create Animation_realistic3 and
       variables_realistic3 folders
All the file associated associated (like .cu and .h files) are in the /src/ folder.

To compile, use the following steps:
1. qrsh -q gpu -l gpu_card=1
2. cd BuddingCode_nucleus5 (note: this is the main folder where you put your budding_cmake folder in, so the name will vary.)
3. cd budding_cmake
4. cd build
5. module list (note: this is to check if you have the correct modules loaded.)
6. module load [...] (note: exact commands skipped because this depends on what modules you need.)
7. make clean; make -j 12;
8. mv virus-model .././Simulation/
9. exit (note: exit the gpu login sesssion when done compiling.)
10. pwd (note: this is just to remind where we are before entering gpu session.)
11. cd ./Simulation/
12. qsub run.sh
(13.) To run interactively, DO NOT exit the gpu session. Instead of exiting, type ./virus-model ../src/Data_structure.xml
