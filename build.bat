g++ -g -Wall -o2 -pthread -static-libstdc++ -std=c++17 -I./include src\differential_corrections_cr3bp.cpp -o DC.exe
:: g++ -g -Wall -o2 -pthread -static-libstdc++ -std=c++17 -I./include -c src/Simulation.cpp -o Simulation.o
g++ -g -Wall -o2 -static-libstdc++ -std=c++17 -I./include -o test.exe src/test.cpp src/Simulation.cpp src/OrbitalElements.cpp src/OrbitalDynamics.cpp 

pause
