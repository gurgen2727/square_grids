CC = gcc
CXX = g++
LINK = g++
MKDIR = mkdir -p
DEL_FILE = rm -f
DEL_DIR = rmdir
SOURCES := $(wildcard src/*.cpp)
OBJECTS := $(patsubst src/%.cpp,obj/%.o,$(SOURCES))
DEPENDS := $(patsubst src/%.cpp,obj/%.dep,$(SOURCES))
bin/program: ./bin $(OBJECTS) bin/out
	echo "./bin/out" > ./bin/program
	echo "octave ./octave/octave.m" >> ./bin/program
	chmod 777 ./bin/program
bin/out: 
	$(LINK) $(OBJECTS) -o bin/out
./bin :
	mkdir -p bin
obj/%.o : src/%.cpp 
	$(CXX) -c $< -o $@
obj/%.dep : src/%.cpp ./obj
	$(CXX) -MM $< -MT "$@ $(putsubst %.dep,%.o,$@)" -MF $@
./obj:
	mkdir -p obj
.PHONY : clean
clean:
	$(DEL_FILE) bin/*
	$(DEL_DIR) bin
	$(DEL_FILE) obj/*
	$(DEL_DIR) obj
-include $(DEPENDS)
