Cube.gif: SI
	convert -delay 5 -loop 0 Sample\ Images/cube* Cube.gif
SI: a.out Test.txt
	./a.out < Test.txt
a.out: cube.cpp
	g++ cube.cpp
	mkdir Sample\ Images
clean:
	rm a.out Cube.gif
	rm -rf Sample\ Images
	