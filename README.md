# triangle
OpenGL triangle visualisation code for triangulation and display of image data.

# build
triangle depends on the following libraries
- [linpng](http://www.libpng.org/pub/png/libpng.html)
- [zlib](http://www.zlib.net/)

I have only tested this code under macos and I link against the Apple provided GLUT (-framework GLUT).

To compile, at the top level directory
```
$ make
```
under macos 10.12 it will throw a bunch of warnings about deprecated OpenGL stuff which I have so far ignored(!)...

The binary will be in `bin/tri`.

Alternatively, I have placed a zip archive in the repo which is a static macos 10.12 build. It *should* work on your system (let me know if it doesn't).

# Usage
triangle is an ongoing work in progress so things change now and then. Right now the app is started on the command line,
```
$ cd bin
$ ./tri [input png image] out.png
```
where `[input png image]` is the name of the input PNG image file you wish to triangulate and `out.png` is the name of a debug image the is generated by the app (it can be called anything you want). This debug image contians the control points generated by the software, displayed as red dots, and the image edges as detected by a Sobol filter. You can ignore this file but a filename *must* be present on the command line.

__NOTE__: currently the input PNG image _must_ be an RGBA png image. I know I should fix this...maybe one day...

To convert a RGB PNG file to an RGBA PNG file I use ImageMagick (or GraphicsMagick). For ImageMagick:
```
$ convert input.png png32:output.png
```
and for GraphicsMagick:
```
$ gm convert input.png png32:output.png
```
should do the trick.

The following options control the visualization

| Key | Action          |
| ------------- | ----------- |
| 1 | white background |
| 2 | black background |
| 3 | greyish background |
| g | grayscale the image (discard color information) |
| = | add contol points (more triangles) |
| - | remove control points (less triangles) |
| [ | increase the percentage of control points the are on an edge|
| ] | decrease the percentage of control points that are on an edge|
| o | decrease the minimum edge value |
| p | increase the minimum edge value |
| r | color triangles greater than some area green |
| , | decrease the area threshold for coloring triangles green  |
| . | increase the area threshold for coloring triangles green |
| SPACE BAR | randomize the control points (you will press space a lot!) |
| f | fade image towards black |
| F | fade image towards white|
| b | toggle border drawing |
| / | toggle triangle drawing mode between GL_TRIANGLES and GL_TRIANGLE_STRIP |
| s | save current image as a png file in the current directory |
| S | save a supersized 300 DPI image in the current directory (takes a little while, but the output is suitable for printing)|
| q | quit |

# Sample output
### Input image
![input sample](https://raw.github.com/freemanjustin/triangle/master/bin/zooey.png)

### sample debug output (out.png in the instructions above) (you can ignore this, and yes it writes it upside down!)
![debug sample2](https://raw.github.com/freemanjustin/triangle/master/bin/debug.png)

### Sample visualisation
![triangle sample3](https://raw.github.com/freemanjustin/triangle/master/bin/sample1.png)
![triangle sample3](https://raw.github.com/freemanjustin/triangle/master/bin/sample2.png)
![triangle sample3](https://raw.github.com/freemanjustin/triangle/master/bin/sample3.png)
![triangle sample3](https://raw.github.com/freemanjustin/triangle/master/bin/sample4.png)
![triangle sample3](https://raw.github.com/freemanjustin/triangle/master/bin/sample5.png)
![triangle sample3](https://raw.github.com/freemanjustin/triangle/master/bin/sample6.png)

# Author
twitter `@_freej_`
instagram `@__freej__`
email freeman.justin@gmail.com
