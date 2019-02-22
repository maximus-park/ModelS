A 3D software package for rasterization, ray-tracing and animation. Build Instructions

Create a directory to build the code:

	$ cd Scotty3D && mkdir build && cd build

Run CMake to generate makefile:

	$ cmake ..

Build the code:

	$ make

Install the executable:

	$ make install

## User Guide for Pencil Sketch

You can run the code on the bunny model by

	$ ./scotty3d bunny.dae

Once you are in MeshEdit mode, press “z” or “Z” on the keyboard to enter Pencil Mode. You will be able to see the pencil rendered bunny now.

You can use the four arrow keys to change the light position, which will create different shadow effects.

You can press “t” or “T” to switch to and back from toon shading mode, which was a progress point for this assignment.

You can press the plus and minus key to change “Cut Off Value” (defined between 1 and 254), which will change the smoothness of bold contour lines. For the bunny model, this value is recommended at 120, which is the default setting. For other models, a value from 90 to 150 will generate a reasonable result. As you increase this value, contour lines tend to become more consistent and contain mostly silhouettes, as more pixels of suggestive contours is being removed from the central area.

The output should look something like this:
Rasterization with brush stroke mode:
![Pencil Sketch](Documentation/images/pencil_sketch.png?raw=true)
Rasterization with toon shading mode:
![Toon Shading](Documentation/images/toon_shade.png?raw=true)

User guide for the entire package is listed under

	Documentation/Scotty3D User Guide _ Computer Graphics _ 15-462_662 Fall 2016.htm
