all:
	make Cube1
	make Cube2
Cube2: a.out input2.txt Images
	./a.out < input2.txt
	convert -delay 5 -loop 0 Images/cube* Cube2.gif
Cube1: a.out input1.txt Images
	./a.out < input1.txt
	convert -delay 5 -loop 0 Images/cube* Cube1.gif
a.out: cube.cpp
	g++ cube.cpp
Images:
	mkdir Images
clean:
	rm a.out *.gif
	rm -rf Images