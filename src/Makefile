CXXFLAGS := -std=c++23 -Wall -Wextra -Wpedantic
LDFLAGS := -l:libdsdp.a -llapack -lblas

all: main visualize tester generate

main:
	g++ $(CXXFLAGS) main.cpp graph.cpp $(LDFLAGS) -o .out/main.out

visualize:
	g++ $(CXXFLAGS) visualize.cpp graph.cpp -lraylib $(LDFLAGS) -o .out/visualize.out

generate:
	g++ $(CXXFLAGS) generate.cpp graph.cpp $(LDFLAGS) -o .out/generate.out

tester:
	g++ $(CXXFLAGS) tester.cpp graph.cpp $(LDFLAGS) -o .out/tester.out