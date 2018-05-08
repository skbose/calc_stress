Solve: main.cpp ./src/SimulatorApp.cpp
	g++ main.cpp ./src/SimulatorApp.cpp -I./includes/ -o Solve -lThea -lboost_system -lboost_thread -lboost_filesystem -lboost_program_options -lm

clean:
	rm Solve
