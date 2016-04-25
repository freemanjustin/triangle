/*

gcc -O3 -g -o tri6 main_with_delauny_sobel.c shader_utils.c tr.c ./delauny_test/Clarkson-Delaunay.c save_image.c -I./delauny_test -framework GLUT -framework OpenGL -lpng -lglew

*/



/**
 * From the OpenGL Programming wikibook: http://en.wikibooks.org/wiki/OpenGL_Programming
 * This file is in the public domain.
 * Contributors: Sylvain Beucler
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
/* Use glew.h instead of gl.h to get all the GL prototypes declared */
#include <GL/glew.h>
/* Using the GLUT library for the base windowing setup */
#include <GLUT/glut.h>
#include "shader_utils.h"

#include <png.h>

#include "tr.h"
#include "loadTexture.h"


//#define TILE_WIDTH 128
//#define TILE_HEIGHT 128
#define TILE_BORDER 10

//#define IMAGE_WIDTH 1024
//#define IMAGE_HEIGHT 768

//#define IMAGE_WIDTH 10240
//#define IMAGE_HEIGHT 7680

//#define FILENAME "supersize.png"

#define INSIDE	1
#define OUTSIDE 0

#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

int TILE_WIDTH;
int TILE_HEIGHT;

int IMAGE_WIDTH;
int IMAGE_HEIGHT;

int width, height;
png_byte color_type;
png_byte bit_depth;

png_structp png_ptr;
png_infop info_ptr;
int number_of_passes;
png_bytep *row_pointers;
png_bytep *out_pointers;


GLuint vbo_triangle, vbo_triangle_colors;
GLuint program;
GLint attribute_coord2d, attribute_v_color, attribute_tex_coord;
GLint uniform_fade;

float cur_fade = 0.5;//sinf(glutGet(GLUT_ELAPSED_TIME) / 1000.0 * (2*3.14) / 5) / 2 + 0.5; // 0->1->0 every 5 seconds
float     edge_precent = 1.0;
int     do_border = 1;
int     do_random = 0;
int     do_wireframe = 0;
int     do_textureMap = 0;
int     do_grayScale = 0;

GLuint  texture[10];


// this is the function to find triangles from a list of points:
unsigned int *BuildTriangleIndexList (void *pointList, float factor, int numberOfInputPoints,
                              int numDimensions, int clockwise, int *numTriangleVertices);

typedef struct{
    float   x;
    float   y;
}pts;

typedef struct{
    float   r;
    float   g;
    float   b;
}vcolor;


typedef struct{
  GLfloat coord2d[2];
  GLfloat v_color[3];
  GLfloat tex_coord[2];
}attributes;

typedef enum{
	ALIAS_MODE_ALIASED,
	ALIAS_MODE_ANTIALIASED,
	ALIAS_MODE_MULTISAMPLE
}AliasMode;


attributes *triangle_attributes;
pts     *polygon;
vcolor  *vertex_color;
int     nTri = 8;
int     nTriVert;
int     nPts = 3000;
int     **edges;
int     min_edge_value = 20;
float   area_threshold = 3000;

AliasMode gMode = ALIAS_MODE_MULTISAMPLE;


void abort_(const char * s, ...)
{
    va_list args;
    va_start(args, s);
    vfprintf(stderr, s, args);
    fprintf(stderr, "\n");
    va_end(args);
    abort();
}

int **malloc2d_int(int dim1, int dim2) {

    int i;
    int **p, *base;

    // Allocate array of pointers and 1D block of integers
    p = (int **)malloc(dim1 * sizeof(int *));
    base = (int *)malloc(dim1 * dim2 * sizeof(int));
    if (p == NULL || base == NULL)
      return NULL;

    // Set pointers for each row
    for (i = 0; i < dim1; i++) {
        p[i] = &base[dim2 * i];
    }

    return p;

}


// r,g,b values are from 0 to 1
// h = [0,360], s = [0,1], v = [0,1]
//		if s == 0, then h = -1 (undefined)
void RGBtoHSV( float r, float g, float b, float *h, float *s, float *v )
{
	float min, max, delta;


	//min = MIN( r, g, b );
	//max = MAX( r, g, b );

    min = r < g ? r : g;
    min = min  < b ? min  : b;

    max = r > g ? r : g;
    max = max  > b ? max  : b;

    *v = max;				// v
	delta = max - min;
	if( max != 0 )
		*s = delta / max;		// s
	else {
		// r = g = b = 0		// s = 0, v is undefined
		*s = 0;
		*h = -1;
		return;
	}
	if( r == max )
		*h = ( g - b ) / delta;		// between yellow & magenta
	else if( g == max )
		*h = 2 + ( b - r ) / delta;	// between cyan & yellow
	else
		*h = 4 + ( r - g ) / delta;	// between magenta & cyan
	*h *= 60;				// degrees
	if( *h < 0 )
		*h += 360;
}
void HSVtoRGB( float *r, float *g, float *b, float h, float s, float v )
{
	int i;
	float f, p, q, t;
	if( s == 0 ) {
		// achromatic (grey)
		*r = *g = *b = v;
		return;
	}
	h /= 60;			// sector 0 to 5
	i = floor( h );
	f = h - i;			// factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );
	switch( i ) {
		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;
		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;
		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;
		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;
		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;
		default:		// case 5:
			*r = v;
			*g = p;
			*b = q;
			break;
	}
}

double r2()
{
    return (double)rand() / (double)RAND_MAX ;
}

double rand_in_range(double min, double max){
  return min + (max - min) * ((double)rand() / (double)RAND_MAX);
}


float GetPixelGray(float r, float g, float b)
{
    return( (( (r) * 0.3 + (g) * 0.59 + (b) * 0.11) ) * 255.0 );
    //return( (( (r) * 0.2126 + (g) * 0.7152 + (b) * 0.0722) ) * 255.0 );

}

void add_noise(){


    int x,y;
    int i,j;
    double rnum;

    png_byte* ptr;
	png_byte* row;

    png_byte* ptrp;
	png_byte* rowp;

    //return;

    //for (i=0; i < dim; ++i) {
    //  ((unsigned char *)surf->pixels)[i] = ((rand() & 1) ? 255 : 0);
    //}

    for (y=0; y<height; y++){
        for (x=0; x<width; x++){

            //printf("y = %d, x = %d, i=%d, j = %d\n", y,x,i,j);
            row = row_pointers[y];
            ptr = &(row[(x)*4]);

            // rnum is either 0 or 255:
            //rnum = ((rand() & 2) ? 255 : 0);

            rnum = rand_in_range(0.0, 1.0);
            //printf("rnum = %f\n", rnum);
            //if(rnum == 255){
            if(rnum > 0.9995){
            printf("rnum = %f\n", rnum);
                ptr[0] = 128.0; //rnum;
                ptr[1] = 128.0; //rnum;
                ptr[2] = 128.0; //rnum;
                ptr[3] = 255;
            }

        }
    }

}

void sobel(int **edges)
{
    int x,y;
    int i,j;

    png_byte* ptr;
	png_byte* row;

    png_byte* ptrp;
	png_byte* rowp;

    float   gray;

    // for sobol filter
    int horizFilter[3][3] = {{ 1,   0,  -1},
                             { 2,   0,  -2},
                             { 1,   0,  -1}};
    int vertFilter[3][3]  = {{ 1,   2,   1},
                             { 0,   0,   0},
                             {-1,  -2,  -1}};
    int pixVal = 0;
    int horizPixVal = 0;
    int vertPixVal = 0;



    if (png_get_color_type(png_ptr, info_ptr) == PNG_COLOR_TYPE_RGB)
        abort_("[process_file] input file is PNG_COLOR_TYPE_RGB but must be PNG_COLOR_TYPE_RGBA "
               "(lacks the alpha channel)");

    if (png_get_color_type(png_ptr, info_ptr) != PNG_COLOR_TYPE_RGBA)
        abort_("[process_file] color_type of input file must be PNG_COLOR_TYPE_RGBA (%d) (is %d)",
               PNG_COLOR_TYPE_RGBA, png_get_color_type(png_ptr, info_ptr));

    // convert to greyscale
    for (y=0; y<height; y++){
        for (x=0; x<width; x++){

            //printf("y = %d, x = %d, i=%d, j = %d\n", y,x,i,j);
            row = row_pointers[y];
            ptr = &(row[(x)*4]);

            gray = GetPixelGray(ptr[0]/255.0,ptr[1]/255.0,ptr[2]/255.0);


            rowp = out_pointers[y];
            ptrp = &(rowp[(x)*4]);

            ptrp[0] = gray;
            ptrp[1] = gray;
            ptrp[2] = gray;
            ptrp[3] = 255;

        }
    }

    // apply sobol filter
    for (y=0; y<height; y++){
        for (x=0; x<width; x++){

            //row = out_pointers[y];
            //ptr = &(row[(x)*4]);

            pixVal = 0;
            horizPixVal = 0;
            vertPixVal = 0;
            if(!((x == 0) || (x == width-1) || (y == 0) || (y == height-1))){    //If the current pixel is along the border, ignore it and set to zero
                for(i = -1; i <= 1; i++){                                               //because the kernel does not align to it
                    for(j = -1; j <= 1; j++){

                        row = out_pointers[y+i];
                        ptr = &(row[(x+j)*4]);

                        horizPixVal += (int)(ptr[0]) * horizFilter[i + 1][j + 1];
                        vertPixVal  += (int)(ptr[0]) * vertFilter[i + 1][j + 1];

                        //horizPixVal += (int)(image[y + j][x + i][0]) * horizFilter[i + 1][j + 1];       //Only need to focus on one of the RGB values since the output is
                        //vertPixVal  += (int)(image[y + j][x + i][0]) * vertFilter[i + 1][j + 1];        //greyscale and all three values are the same
                    }
                }
            }
            pixVal = sqrt((horizPixVal * horizPixVal) + (vertPixVal * vertPixVal));     //Calculate magnitude
            pixVal = sqrt(horizPixVal * horizPixVal);
            if(pixVal > 255) pixVal = 255;                                              //Clamp value within 8-bit range
            //filteredImage[y][x][0] = (unsigned char)pixVal;

            //if(pixVal < 150) pixVal = 0;


            //if(pixVal > 250) pixVal = 255;

            //row = out_pointers[y];
            //ptr = &(row[(x)*4]);
            //printf("pixVal = %d, gray = %d\n", pixVal, ptr[0]);

            /*
            row = row_pointers[y];
            ptr = &(row[(x)*4]);

            ptr[0] = pixVal;
            ptr[1] = pixVal;
            ptr[2] = pixVal;
            ptr[3] = 255.0;
            */
            edges[(height-1)-y][x] = pixVal;

        }
    }


    // copy grayscale edge image into the output image memory locations (rowp and ptrp)
    for (y=0; y<height; y++){
        for (x=0; x<width; x++){
            // output image pointers
            //row = row_pointers[y];
            //ptr = &(row[x*4]);

            // processed image pointers
            rowp = out_pointers[y];
            ptrp = &(rowp[(x)*4]);

            ptrp[0] = edges[y][x];//ptr[0];
            ptrp[1] = edges[y][x];//ptr[1];
            ptrp[2] = edges[y][x];//ptr[2];
            ptrp[3] = 255.0;//ptr[3];

            //printf("Output at position [ %d - %d ] has RGBA values: %d - %d - %d - %d. Magnitude is: %f\n",
            //       x, y, ptr[0], ptr[1], ptr[2], ptr[3], (sqrt(ptr[0]/255.0*ptr[0]/255.0+ptr[1]/255.0*ptr[1]/255.0+ptr[2]/255.0*ptr[2]/255.0)));
        }
    }
}




int process_triangle_poly(pts *polygon, vcolor *vertex_color){

    int i;
    int x,y;
    float xmin,xmax,ymin,ymax;
    int nCoords;
    int count;
    float   r,g,b;
    float   h_col,s_col,v_col;
    float   gray;
    float   ave_r, ave_g, ave_b;

    float   len_a, len_b, len_c, s, area;

    pts point;
    png_byte* ptr;
	png_byte* row;

    nCoords = 3;
    count = 0;
    r = g = b = 0.0;

    xmax = -1;
    xmin = 9999999;
    ymax = -1;
    ymin = 9999999;


    for(i=0;i<nCoords;i++){
        if(polygon[i].x > (float)xmax)
            xmax = (int)polygon[i].x;
        if(polygon[i].x < (float)xmin)
            xmin = (int)polygon[i].x;
        if(polygon[i].y > (float)ymax)
            ymax = (int)(polygon[i].y);
        if(polygon[i].y < (float)ymin)
            ymin = (int)(polygon[i].y);
    }

    //for(y=0;y<height;y++){
    //    for(x=0;x<width;x++){
    for(y=height-ymax;y<height-ymin;y++){
        for(x=xmin;x<xmax;x++){

            point.x = (float)x;
            point.y = height-(float)y;
            if(pnpoly(polygon, nCoords, point)){    // if true then pixel is inside the triangle

                //printf("%d %d is INSIDE\n", i,j);
                row = row_pointers[y];
                ptr = &(row[(x)*4]);

                if(do_grayScale == 1){

                    gray = (int)(0.2126 * (float)ptr[0]+ 0.7152 * (float)ptr[1] + 0.0722 * (float)ptr[2]);
                    r += gray;
                    g += gray;
                    b += gray;
                    //printf("gray it: %f, %f, %f\n", vertex_color->r, vertex_color->g, vertex_color->b);
                }
                else{
                    r += (float)ptr[0];
                    g += (float)ptr[1];
                    b += (float)ptr[2];
                }
                count++;
            }
        }
    }

    if(count == 0){
        //printf("got nothing\n");
        return -1;
    }

    if(do_random == 1){

        // only modify the color if the triangle area is above a threshold

        // triangle area:

        //      s = (a+b+c)/2;
        //      area = sqrt(s*(s-a)*(s-b)*(s-c));
        //
        //  where a, b and c are the lengths of the sides of the triangle
        //
        //  Length of each side is:

        //      distance between x(x1,y1) and y(x2,y2) is:

        //      d = sqrt( (x2 - x1) + (y2 - y1) )

        //  first get line lengths


        len_a = sqrt (pow((polygon[1].x -polygon[0].x),2.0) + pow((polygon[1].y - polygon[0].y),2.0) );
        len_b = sqrt (pow((polygon[2].x -polygon[1].x),2.0) + pow((polygon[2].y - polygon[1].y),2.0) );
        len_c = sqrt (pow((polygon[2].x -polygon[0].x),2.0) + pow((polygon[2].y - polygon[0].y),2.0) );

        /*
        printf("a: (%f, %f) to (%f, %f): distane is %f\n", polygon[0].x, polygon[0].y, polygon[1].x, polygon[1].y, len_a);

        printf("b: (%f, %f) to (%f, %f): distane is %f\n", polygon[1].x, polygon[1].y, polygon[2].x, polygon[2].y, len_b);

        printf("c: (%f, %f) to (%f, %f): distane is %f\n", polygon[2].x, polygon[2].y, polygon[0].x, polygon[0].y, len_c);
        */

        s = (len_a + len_b + len_c)/2.0;

        area = sqrt(s*(s-len_a)*(s-len_b)*(s-len_c));

        //printf("area = %f\n\n",area);

        if(area > area_threshold){
        //if(rand_in_range(0.0,1.0)>0.90){

            //printf("doing test: area = %f\n", area);

            //if(rand_in_range(0.0,1.0)>0.5){
            {

                // for n_woods
                // this will give a nice orange color
                /*
                vertex_color->r = rand_in_range(0.70,0.98);
                vertex_color->g = 0.4;
                vertex_color->b = 0.1;
                */

                // for n_woods
                // blue color
                //vertex_color->r = 40.0/255.0;//rand_in_range(0.70,0.98);
                //vertex_color->g = 83.0/255.0;
                //vertex_color->b = rand_in_range(80.0/255.0,110.0/255.0);//255.0/255.0;155

                // nwoods greens
                /*
                vertex_color->r = 103/255.0;
                vertex_color->g = rand_in_range(180,235)/255.0;//205/255.0;
                vertex_color->b = 35/255.0;
                */


                // this is the equivalent green to the above
                // experimental rgb -> hsv -> randomize over value -> rgb
                RGBtoHSV( 103.0/255.0, 205.0/255.0, 35.0/255.0, &h_col, &s_col, &v_col );
                HSVtoRGB(&vertex_color->r, &vertex_color->g, &vertex_color->b, h_col, s_col, v_col+rand_in_range(-0.3,0.1));

                // nwoods blue
                //RGBtoHSV( 75/255.0, 94/255.0, 214/255.0, &h_col, &s_col, &v_col );
                //HSVtoRGB(&vertex_color->r, &vertex_color->g, &vertex_color->b, h_col, s_col, v_col+rand_in_range(-0.3,0.1));

                // nwoods red
                //RGBtoHSV( 237/255.0, 60/255.0, 47/255.0, &h_col, &s_col, &v_col );
                //HSVtoRGB(&vertex_color->r, &vertex_color->g, &vertex_color->b, h_col, s_col, v_col+rand_in_range(-0.2,0.2));


                // nwoods orange
                //RGBtoHSV( 184/255.0, 82/255.0, 20/255.0, &h_col, &s_col, &v_col );
                //HSVtoRGB(&vertex_color->r, &vertex_color->g, &vertex_color->b, h_col, s_col, v_col+rand_in_range(-0.3,0.3));


                // nwoods yellows
                /*
                vertex_color->r = rand_in_range(220,245)/255.0;
                vertex_color->g = rand_in_range(220,255)/255.0;//205/255.0;
                vertex_color->b = 0.0/255.0;//rand_in_range(1,5)/255.0;
                */


                // for grey dino
                // red color
                /*
                vertex_color->r = rand_in_range(0.8,0.98);
                vertex_color->g = 0.1;
                vertex_color->b = 0.05;
                */
            }
        }
        else{
            //printf("I found %d pixels inside this triangle\n", count);
            // average the pixel values
            vertex_color->r = (r/(float)(count))/255.0;
            vertex_color->g = (g/(float)(count))/255.0;
            vertex_color->b = (b/(float)(count))/255.0;
            //printf("colors are: %f, %f, %f\n", vertex_color->r, vertex_color->g, vertex_color->b);


            if(do_grayScale == 1){
                RGBtoHSV( vertex_color->r, vertex_color->g, vertex_color->b, &h_col, &s_col, &v_col );
                HSVtoRGB(&vertex_color->r, &vertex_color->g, &vertex_color->b, h_col, s_col, v_col+rand_in_range(-0.2,0.0));
            }
        }
    }
    else if(do_wireframe == 1){
        // if wireframe is true the set the color to black
        vertex_color->r = 0.0/255.0;
        vertex_color->g = 0.0/255.0;
        vertex_color->b = 0.0/255.0;
    }
    else{
        vertex_color->r = (r/(float)(count))/255.0;
        vertex_color->g = (g/(float)(count))/255.0;
        vertex_color->b = (b/(float)(count))/255.0;

        if(do_grayScale == 1){
            RGBtoHSV( vertex_color->r, vertex_color->g, vertex_color->b, &h_col, &s_col, &v_col );
            HSVtoRGB(&vertex_color->r, &vertex_color->g, &vertex_color->b, h_col, s_col, v_col+rand_in_range(-0.1,0.1));
        }
    }

    return 1;

}


int pnpoly(pts *polygon, int npol, pts p){

	int i;
	int c = 0;

	int	blocksize = 8;
	int	blocklimit;

	blocklimit = (npol / blocksize) * blocksize;

	// special wrap around case:
	if ((((polygon[0].y <= p.y) && (p.y < polygon[npol-1].y )) ||
		(( polygon[npol-1].y <= p.y) && (p.y < polygon[0].y ))) &&
		(p.x < (polygon[npol-1].x - polygon[0].x) * (p.y - polygon[0].y) / (polygon[npol-1].y - polygon[0].y) + polygon[0].x))
   		       		c = !c;

	i = 1;
	while(i<blocklimit){
		if ((((polygon[i].y <= p.y) && (p.y < polygon[i-1].y )) ||
				(( polygon[i-1].y <= p.y) && (p.y < polygon[i].y ))) &&
				(p.x < (polygon[i-1].x - polygon[i].x) * (p.y - polygon[i].y) / (polygon[i-1].y - polygon[i].y) + polygon[i].x))
   		       		c = !c;
		if ((((polygon[i+1].y <= p.y) && (p.y < polygon[i].y )) ||
				(( polygon[i].y <= p.y) && (p.y < polygon[i+1].y ))) &&
				(p.x < (polygon[i].x - polygon[i+1].x) * (p.y - polygon[i+1].y) / (polygon[i].y - polygon[i+1].y) + polygon[i+1].x))
   		       		c = !c;
		if ((((polygon[i+2].y <= p.y) && (p.y < polygon[i+1].y )) ||
				(( polygon[i+1].y <= p.y) && (p.y < polygon[i+2].y ))) &&
				(p.x < (polygon[i+1].x - polygon[i+2].x) * (p.y - polygon[i+2].y) / (polygon[i+1].y - polygon[i+2].y) + polygon[i+2].x))
   		       		c = !c;
		if ((((polygon[i+3].y <= p.y) && (p.y < polygon[i+2].y )) ||
				(( polygon[i+2].y <= p.y) && (p.y < polygon[i+3].y ))) &&
				(p.x < (polygon[i+2].x - polygon[i+3].x) * (p.y - polygon[i+3].y) / (polygon[i+2].y - polygon[i+3].y) + polygon[i+3].x))
   		       		c = !c;
		if ((((polygon[i+4].y <= p.y) && (p.y < polygon[i+3].y )) ||
				(( polygon[i+3].y <= p.y) && (p.y < polygon[i+4].y ))) &&
				(p.x < (polygon[i+3].x - polygon[i+4].x) * (p.y - polygon[i+4].y) / (polygon[i+3].y - polygon[i+4].y) + polygon[i+4].x))
   		       		c = !c;
		if ((((polygon[i+5].y <= p.y) && (p.y < polygon[i+4].y )) ||
				(( polygon[i+4].y <= p.y) && (p.y < polygon[i+5].y ))) &&
				(p.x < (polygon[i+4].x - polygon[i+5].x) * (p.y - polygon[i+5].y) / (polygon[i+4].y - polygon[i+5].y) + polygon[i+5].x))
   		       		c = !c;
		if ((((polygon[i+6].y <= p.y) && (p.y < polygon[i+5].y )) ||
				(( polygon[i+5].y <= p.y) && (p.y < polygon[i+6].y ))) &&
				(p.x < (polygon[i+5].x - polygon[i+6].x) * (p.y - polygon[i+6].y) / (polygon[i+5].y - polygon[i+6].y) + polygon[i+6].x))
   		       		c = !c;
		if ((((polygon[i+7].y <= p.y) && (p.y < polygon[i+6].y )) ||
				(( polygon[i+6].y <= p.y) && (p.y < polygon[i+7].y ))) &&
				(p.x < (polygon[i+6].x - polygon[i+7].x) * (p.y - polygon[i+7].y) / (polygon[i+6].y - polygon[i+7].y) + polygon[i+7].x))
   		       		c = !c;

		i+=blocksize;
	}

	// mop up the remainder...
	if(i<npol){
		do{
			if ((((polygon[i].y <= p.y) && (p.y < polygon[i-1].y )) ||
				(( polygon[i-1].y <= p.y) && (p.y < polygon[i].y ))) &&
				(p.x < (polygon[i-1].x - polygon[i].x) * (p.y - polygon[i].y) / (polygon[i-1].y - polygon[i].y) + polygon[i].x))
   		       		c = !c;
			i++;
		}while(i<npol);
	}

	return c;
}




void read_png_file(char* file_name)
{
    int y;
    char header[8];    // 8 is the maximum size that can be checked

    /* open file and test for it being a png */
    FILE *fp = fopen(file_name, "rb");
    if (!fp)
        abort_("[read_png_file] File %s could not be opened for reading", file_name);
    fread(header, 1, 8, fp);
    if (png_sig_cmp(header, 0, 8))
        abort_("[read_png_file] File %s is not recognized as a PNG file", file_name);


    /* initialize stuff */
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
        abort_("[read_png_file] png_create_read_struct failed");

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
        abort_("[read_png_file] png_create_info_struct failed");

    if (setjmp(png_jmpbuf(png_ptr)))
        abort_("[read_png_file] Error during init_io");

    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);

    png_read_info(png_ptr, info_ptr);

    width = png_get_image_width(png_ptr, info_ptr);
    height = png_get_image_height(png_ptr, info_ptr);
    color_type = png_get_color_type(png_ptr, info_ptr);
    bit_depth = png_get_bit_depth(png_ptr, info_ptr);

    number_of_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);


    /* read file */
    if (setjmp(png_jmpbuf(png_ptr)))
        abort_("[read_png_file] Error during read_image");

    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (y=0; y<height; y++)
        row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));

    png_read_image(png_ptr, row_pointers);

    out_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (y=0; y<height; y++)
        out_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));

    //png_read_image(png_ptr, out_pointers);

    fclose(fp);
}


void write_png_file(char* file_name)
{
    int y;

    /* create file */
    FILE *fp = fopen(file_name, "wb");
    if (!fp)
        abort_("[write_png_file] File %s could not be opened for writing", file_name);


    /* initialize stuff */
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
        abort_("[write_png_file] png_create_write_struct failed");

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
        abort_("[write_png_file] png_create_info_struct failed");

    if (setjmp(png_jmpbuf(png_ptr)))
        abort_("[write_png_file] Error during init_io");

    png_init_io(png_ptr, fp);


    /* write header */
    if (setjmp(png_jmpbuf(png_ptr)))
        abort_("[write_png_file] Error during writing header");

    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);


    /* write bytes */
    if (setjmp(png_jmpbuf(png_ptr)))
        abort_("[write_png_file] Error during writing bytes");

    png_write_image(png_ptr, out_pointers);


    /* end write */
    if (setjmp(png_jmpbuf(png_ptr)))
        abort_("[write_png_file] Error during end of write");

    png_write_end(png_ptr, NULL);

    /* cleanup heap allocation */
    for (y=0; y<height; y++)
        free(row_pointers[y]);
    free(row_pointers);

    for (y=0; y<height; y++)
        free(out_pointers[y]);
    free(out_pointers);

    fclose(fp);
}


int init_resources()
{

  glGenBuffers(1, &vbo_triangle);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle);
  //glBufferData(GL_ARRAY_BUFFER, sizeof(triangle_attributes), triangle_attributes, GL_STATIC_DRAW);
  glBufferData(GL_ARRAY_BUFFER, sizeof(attributes)*nTriVert, triangle_attributes, GL_STATIC_DRAW);

  GLint link_ok = GL_FALSE;

  GLuint vs, fs;
  if ((vs = create_shader("shader/triangle.v.glsl", GL_VERTEX_SHADER))   == 0) return 0;
  if ((fs = create_shader("shader/triangle.f.glsl", GL_FRAGMENT_SHADER)) == 0) return 0;

  program = glCreateProgram();
  glAttachShader(program, vs);
  glAttachShader(program, fs);
  glLinkProgram(program);
  glGetProgramiv(program, GL_LINK_STATUS, &link_ok);
  if (!link_ok) {
    fprintf(stderr, "glLinkProgram:");
    print_log(program);
    return 0;
  }

  const char* attribute_name;
  attribute_name = "coord2d";
  attribute_coord2d = glGetAttribLocation(program, attribute_name);
  if (attribute_coord2d == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return 0;
  }
  attribute_name = "v_color";
  attribute_v_color = glGetAttribLocation(program, attribute_name);
  if (attribute_v_color == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return 0;
  }


  /*
  attribute_name = "tex_coord";
  attribute_tex_coord = glGetAttribLocation(program, attribute_name);
  if (attribute_tex_coord == -1) {
    fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
    return 0;
  }
  */

  const char* uniform_name;
  uniform_name = "fade";
  uniform_fade = glGetUniformLocation(program, uniform_name);
  if (uniform_fade == -1) {
    fprintf(stderr, "Could not bind uniform_fade %s\n", uniform_name);
    return 0;
  }

  return 1;
}

void update_resources_wireframe(){

    int i;

    for(i=0;i<nTriVert;i+=3){   //for each triangle

        // rescale triangle coords to match image dimensions
        // points are in -1 <= x <= 1,     -1 <= y <= 1
        // image is       0 <= x <= width,  0 <= y <= height

        /*
        polygon[0].x = width * (0.5*triangle_attributes[i].coord2d[0] + 0.5);
        polygon[0].y = height *(0.5*triangle_attributes[i].coord2d[1] + 0.5);

        polygon[1].x = width * (0.5*triangle_attributes[i+1].coord2d[0] + 0.5);
        polygon[1].y = height *(0.5*triangle_attributes[i+1].coord2d[1] + 0.5);

        polygon[2].x = width * (0.5*triangle_attributes[i+2].coord2d[0] + 0.5);
        polygon[2].y = height * (0.5*triangle_attributes[i+2].coord2d[1] + 0.5);
        */

        /*
        printf("polygon in image space:\n");
        printf("tri = %d\n", i/3);
        printf("\tvertex: %f, %f (was %f, %f)\n", polygon[0].x, polygon[0].y, 0.5*triangle_attributes[i].coord2d[0] + 0.5, 0.5*triangle_attributes[i].coord2d[1] + 0.5);
        printf("\tvertex: %f, %f (was %f, %f)\n", polygon[1].x, polygon[1].y, 0.5*triangle_attributes[i+1].coord2d[0]+0.5,0.5*triangle_attributes[i+1].coord2d[1]+0.5);
        printf("\tvertex: %f, %f (was %f, %f)\n", polygon[2].x, polygon[2].y, 0.5*triangle_attributes[i+2].coord2d[0]+0.5, 0.5*triangle_attributes[i+2].coord2d[1]+0.5);
        */

        process_triangle_poly(polygon, vertex_color);    // find out what pixels are in this triangle

        triangle_attributes[i].v_color[0] = vertex_color->r;
        triangle_attributes[i].v_color[1] = vertex_color->g;
        triangle_attributes[i].v_color[2] = vertex_color->b;

        triangle_attributes[i+1].v_color[0] = vertex_color->r;
        triangle_attributes[i+1].v_color[1] = vertex_color->g;
        triangle_attributes[i+1].v_color[2] = vertex_color->b;

        triangle_attributes[i+2].v_color[0] = vertex_color->r;
        triangle_attributes[i+2].v_color[1] = vertex_color->g;
        triangle_attributes[i+2].v_color[2] = vertex_color->b;
    }

    glBufferData(GL_ARRAY_BUFFER, nTriVert*sizeof(attributes), triangle_attributes, GL_STATIC_DRAW);


}


void update_resources(){
    int i,l,x,y;
    double  rnum;
    double  r1,r2;
    float     *points;
    float   fac = RAND_MAX;
    int     nDim = 2;
    int     cwise = 1;
    pts     *point_list;


    unsigned int    *triIndexList;
    static int      done_load_texture = 0;
    static          pngInfo infoLayer;

    points = malloc(nDim*nPts*sizeof(float));
    point_list = malloc(nPts * sizeof(pts));


    if(do_textureMap == 1){

        // load textures only once
        if(done_load_texture == 0){

            /*
            // open the texture map image
            texture[0] = pngBind("texture_test.png",
									PNG_BUILDMIPMAPS,
									PNG_SOLID,
									&infoLayer,
									GL_CLAMP_TO_EDGE,
									GL_LINEAR_MIPMAP_LINEAR,
									GL_LINEAR_MIPMAP_LINEAR);

            printf("done load texture\n");
            */
            done_load_texture = 1;
        }

    }

    l = 0;

    if(do_border == 1){
        for(i=0;i<=10;i+=2){
            points[l] = (float)i/10.0;
            if( (i!=0) && (i!=10) )
                points[l] += rand_in_range(-0.1,0.1);
            l++;
            points[l] = 0.0; l++;
        }

        for(i=1;i<=10;i+=2){
            points[l] = (float)i/10.0;
            if( (i!=0) && (i!=10) )
                points[l] += rand_in_range(-0.1,0.1);
            l++;
            points[l] = 1.0; l++;
        }

        for(i=0;i<=10;i+=2){
            points[l] = 0.0; l++;
            points[l] = (float)i/10.0;
            if( (i!=0) && (i!=10) )
                points[l] += rand_in_range(-0.1,0.1);
            l++;
        }

        for(i=0;i<=10;i+=2){
            points[l] = 1.0; l++;
            points[l] = (float)i/10.0;
            if( (i!=0) && (i!=10) )
                points[l] += rand_in_range(-0.1,0.1);
            l++;
        }
    }

    // how may pixels are in the image that represent an edge?
    //l = 0;
    do{
        r1 = rand_in_range(0.0,1.0);
        r2 = rand_in_range(0.0,1.0);
        // convert random numbers to image coordinates
        y = (int)((float)height * r1);
        x = (int)((float)width  * r2);
        // is this on an adge?
        if(edges[y][x] > min_edge_value){   // found an edge
            points[l] = r2; l++;
            points[l] = r1; l++;
        }

    }while(l<(nDim*nPts)*edge_precent);

    for(i=l;i<2*nPts;i++){
        points[i] = rand_in_range(0.0,1.0);
    }

    /*
    for(i=0;i<=10;i+=2){
        points[l] = (float)i/10.0; l++;
        points[l] = 0.0; l++;
    }

    for(i=1;i<=10;i+=2){
        points[l] = (float)i/10.0; l++;
        points[l] = 1.0; l++;
    }

    for(i=0;i<=10;i+=2){
        points[l] = 0.0; l++;
        points[l] = (float)i/10.0; l++;
    }

    for(i=0;i<=10;i+=2){
        points[l] = 1.0; l++;
        points[l] = (float)i/10.0; l++;
    }
    */

    for(i=l;i<2*nPts;i++){
        points[i] = rand_in_range(0.0,1.0);
    }



    l = 0;
    for(i=0;i<2*nPts;i+=2){
        point_list[l].x = points[i];
        point_list[l].y = points[i+1];
        l++;
    }

    //printf("#building triangle list\n");
    triIndexList = BuildTriangleIndexList((void*)points, fac, nPts, nDim, cwise, &nTriVert);


    nTri = (nTriVert/3);
    free(triangle_attributes);
    triangle_attributes = malloc(nTriVert * sizeof(attributes));

    //printf("#done\n");

    //printf("#numTriangle Vert = %d\n", nTriVert);
    //printf("#nTriangles = %d\n", nTriVert/3);


    /*
    l = 0;
    for(i=0;i<2*nPts;i+=2){
        printf("#point[ %d ]= %f,%f\n", l, points[i], points[i+1]);
        l++;
    }
    */

    /*
    l = 0;
    for(i=0;i<nPts;i++){
        printf("#point_list[ %d ] = %f,%f\n", i, point_list[i].x, point_list[i].y);
        l++;
    }
    */

    // set up triangles for display
    //nTri = (nTriVert/3);

    for(i=0;i<nTriVert;i+=3){
        //printf("#triangle %d hs vertices = %d, %d, %d\n",l, triIndexList[i],triIndexList[i+1],triIndexList[i+2]);
        /*
        r1 = rand_in_range(0.0, 0.1);
        r2 = rand_in_range(0.4, 1.0);
        rnum = rand_in_range(0.4, 1.0);


        triangle_attributes[i].v_color[0] = rnum; //r1; //0.31640625 + r1;
        triangle_attributes[i].v_color[1] = rnum ; //r2; //0.65625 + r2;
        triangle_attributes[i].v_color[2] = rnum;

        triangle_attributes[i+1].v_color[0] = rnum;//r1; //0.31640625 + r1;
        triangle_attributes[i+1].v_color[1] = rnum;//r2; //0.65625 + r2;
        triangle_attributes[i+1].v_color[2] = rnum;

        triangle_attributes[i+2].v_color[0] = rnum;//r1; //0.31640625 + r1;
        triangle_attributes[i+2].v_color[1] = rnum;//r2; //0.65625 + r2;
        triangle_attributes[i+2].v_color[2] = rnum;
        */

        triangle_attributes[i].coord2d[0] = (point_list[triIndexList[i]].x - 0.5)/0.5;
        triangle_attributes[i].coord2d[1] = (point_list[triIndexList[i]].y - 0.5)/0.5;

        triangle_attributes[i+1].coord2d[0] = (point_list[triIndexList[i+1]].x - 0.5)/0.5;
        triangle_attributes[i+1].coord2d[1] = (point_list[triIndexList[i+1]].y - 0.5)/0.5;

        triangle_attributes[i+2].coord2d[0] = (point_list[triIndexList[i+2]].x - 0.5)/0.5;
        triangle_attributes[i+2].coord2d[1] = (point_list[triIndexList[i+2]].y - 0.5)/0.5;


    }

    for(i=0;i<nTriVert;i+=3){   //for each triangle

        // rescale triangle coords to match image dimensions
        // points are in -1 <= x <= 1,     -1 <= y <= 1
        // image is       0 <= x <= width,  0 <= y <= height


        polygon[0].x = width * (0.5*triangle_attributes[i].coord2d[0] + 0.5);
        polygon[0].y = height *(0.5*triangle_attributes[i].coord2d[1] + 0.5);

        polygon[1].x = width * (0.5*triangle_attributes[i+1].coord2d[0] + 0.5);
        polygon[1].y = height *(0.5*triangle_attributes[i+1].coord2d[1] + 0.5);

        polygon[2].x = width * (0.5*triangle_attributes[i+2].coord2d[0] + 0.5);
        polygon[2].y = height * (0.5*triangle_attributes[i+2].coord2d[1] + 0.5);

        /*
        printf("polygon in image space:\n");
        printf("tri = %d\n", i/3);
        printf("\tvertex: %f, %f (was %f, %f)\n", polygon[0].x, polygon[0].y, 0.5*triangle_attributes[i].coord2d[0] + 0.5, 0.5*triangle_attributes[i].coord2d[1] + 0.5);
        printf("\tvertex: %f, %f (was %f, %f)\n", polygon[1].x, polygon[1].y, 0.5*triangle_attributes[i+1].coord2d[0]+0.5,0.5*triangle_attributes[i+1].coord2d[1]+0.5);
        printf("\tvertex: %f, %f (was %f, %f)\n", polygon[2].x, polygon[2].y, 0.5*triangle_attributes[i+2].coord2d[0]+0.5, 0.5*triangle_attributes[i+2].coord2d[1]+0.5);
        */

        process_triangle_poly(polygon, vertex_color);    // find out what pixels are in this triangle

        triangle_attributes[i].v_color[0] = vertex_color->r;
        triangle_attributes[i].v_color[1] = vertex_color->g;
        triangle_attributes[i].v_color[2] = vertex_color->b;

        triangle_attributes[i+1].v_color[0] = vertex_color->r;
        triangle_attributes[i+1].v_color[1] = vertex_color->g;
        triangle_attributes[i+1].v_color[2] = vertex_color->b;

        triangle_attributes[i+2].v_color[0] = vertex_color->r;
        triangle_attributes[i+2].v_color[1] = vertex_color->g;
        triangle_attributes[i+2].v_color[2] = vertex_color->b;

        // texture coordinates
        triangle_attributes[i].tex_coord[0] = 0.5;
        triangle_attributes[i].tex_coord[1] = 1.0;
        triangle_attributes[i+1].tex_coord[0] = 0.0;
        triangle_attributes[i+1].tex_coord[1] = 0.0;
        triangle_attributes[i+2].tex_coord[0] = 1.0;
        triangle_attributes[i+2].tex_coord[1] = 0.0;
    }

    glBufferData(GL_ARRAY_BUFFER, nTriVert*sizeof(attributes), triangle_attributes, GL_STATIC_DRAW);

}

void onIdle() {
  //float cur_fade = 0.8;//sinf(glutGet(GLUT_ELAPSED_TIME) / 1000.0 * (2*3.14) / 5) / 2 + 0.5; // 0->1->0 every 5 seconds
  //glUseProgram(program);
  //glUniform1f(uniform_fade, cur_fade);
  //glutPostRedisplay();
}

void onDisplay()
{
  int i;
  //glClearColor(0.0, 0.0, 0.0, 1.0);
  //glClear(GL_COLOR_BUFFER_BIT);

  glClear(GL_COLOR_BUFFER_BIT);
  glClear(GL_DEPTH_BUFFER_BIT);

  glUseProgram(program);
  glUniform1f(uniform_fade, cur_fade);

  /*
  //Start alias mode
	switch( gMode )
	{
	    case ALIAS_MODE_ALIASED:
            glDisable( GL_LINE_SMOOTH );
			glDisable( GL_POLYGON_SMOOTH );
			glDisable( GL_MULTISAMPLE );
            break;

		case ALIAS_MODE_ANTIALIASED:
			glEnable( GL_LINE_SMOOTH );
			glEnable( GL_POLYGON_SMOOTH );
            glDisable( GL_MULTISAMPLE );
			break;

		case ALIAS_MODE_MULTISAMPLE:
            glDisable( GL_LINE_SMOOTH );
			glDisable( GL_POLYGON_SMOOTH );
			glEnable( GL_MULTISAMPLE );
            glEnable( GL_MULTISAMPLE_ARB );
			break;
	}
    */

  //glUseProgram(program);

  glEnableVertexAttribArray(attribute_coord2d);
  glEnableVertexAttribArray(attribute_v_color);

  /*
  if(do_textureMap == 1){
  glEnableVertexAttribArray(attribute_tex_coord);
  }
  */


  glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle);
  glVertexAttribPointer(
    attribute_coord2d,   // attribute
    2,                   // number of elements per vertex, here (x,y)
    GL_FLOAT,            // the type of each element
    GL_FALSE,            // take our values as-is
    sizeof(attributes),  // next coord2d appears every 5 floats
    0                    // offset of first element
  );
  glVertexAttribPointer(
    attribute_v_color,      // attribute
    3,                      // number of elements per vertex, here (r,g,b)
    GL_FLOAT,               // the type of each element
    GL_FALSE,               // take our values as-is
    sizeof(attributes),  // stride
    //(GLvoid*) (2 * sizeof(GLfloat))     // offset of first element
    (GLvoid*) offsetof(attributes, v_color)  // offset
  );


/*
if(do_textureMap == 1){
  glVertexAttribPointer(
    attribute_tex_coord,      // attribute
    2,                      // number of elements per vertex, here (r,g,b)
    GL_FLOAT,               // the type of each element
    GL_FALSE,               // take our values as-is
    sizeof(attributes),  // stride
    //(GLvoid*) (2 * sizeof(GLfloat))     // offset of first element
    (GLvoid*) offsetof(attributes, tex_coord)  // offset
  );
}
*/


/*
  // Push each element in buffer_vertices to the vertex shader
  glDrawArrays(GL_TRIANGLES, 0, 3*nTri);

// Turn on wireframe mode
glPolygonMode(GL_FRONT, GL_LINE);
glPolygonMode(GL_BACK, GL_LINE);

// Draw the tri
glDrawArrays(GL_TRIANGLES, 0, 3*nTri);

// Turn off wireframe mode
glPolygonMode(GL_FRONT, GL_FILL);
glPolygonMode(GL_BACK, GL_FILL);
*/

    //case ALIAS_MODE_MULTISAMPLE:
            glDisable( GL_LINE_SMOOTH );
			glDisable( GL_POLYGON_SMOOTH );
			glEnable( GL_MULTISAMPLE );
            glEnable( GL_MULTISAMPLE_ARB );

    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glDrawArrays(GL_TRIANGLES, 0, 3*nTri);
    //glDrawArrays(GL_TRIANGLE_STRIP,0,3*nTri);


    //case ALIAS_MODE_ANTIALIASED:
			glEnable( GL_LINE_SMOOTH );
			glEnable( GL_POLYGON_SMOOTH );
            glDisable( GL_MULTISAMPLE );

    glLineWidth(1.0f);
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    glDrawArrays(GL_TRIANGLES, 0, 3*nTri);
    //glDrawArrays(GL_TRIANGLE_STRIP,0,3*nTri);

  glDisableVertexAttribArray(attribute_coord2d);
  glDisableVertexAttribArray(attribute_v_color);
  /*
  if(do_textureMap == 1){
    glDisableVertexAttribArray(attribute_tex_coord);
  }
  */

  glutSwapBuffers();
}



void draw_tr()
{
  int i;

  glClear(GL_COLOR_BUFFER_BIT);
  glClear(GL_DEPTH_BUFFER_BIT);


  glEnableVertexAttribArray(attribute_coord2d);
  glEnableVertexAttribArray(attribute_v_color);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_triangle);
  glVertexAttribPointer(
    attribute_coord2d,   // attribute
    2,                   // number of elements per vertex, here (x,y)
    GL_FLOAT,            // the type of each element
    GL_FALSE,            // take our values as-is
    sizeof(attributes),  // next coord2d appears every 5 floats
    0                    // offset of first element
  );
  glVertexAttribPointer(
    attribute_v_color,      // attribute
    3,                      // number of elements per vertex, here (r,g,b)
    GL_FLOAT,               // the type of each element
    GL_FALSE,               // take our values as-is
    sizeof(attributes),  // stride
    //(GLvoid*) (2 * sizeof(GLfloat))     // offset of first element
    (GLvoid*) offsetof(attributes, v_color)  // offset
  );

    //case ALIAS_MODE_MULTISAMPLE:
            glDisable( GL_LINE_SMOOTH );
			glDisable( GL_POLYGON_SMOOTH );
			glEnable( GL_MULTISAMPLE );
            glEnable( GL_MULTISAMPLE_ARB );

    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glDrawArrays(GL_TRIANGLES, 0, 3*nTri);

    //case ALIAS_MODE_ANTIALIASED:
			glEnable( GL_LINE_SMOOTH );
			glEnable( GL_POLYGON_SMOOTH );
            glDisable( GL_MULTISAMPLE );

    //glLineWidth(1.0f);
    // for supersized images make the lines thicker
    glLineWidth(1.0f);
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    glDrawArrays(GL_TRIANGLES, 0, 3*nTri);

  glDisableVertexAttribArray(attribute_coord2d);
  glDisableVertexAttribArray(attribute_v_color);


}



/* Do a demonstration of tiled rendering */
static void tile_Display(void)
{
   TRcontext *tr;
   GLubyte *buffer;
   GLubyte *tile;
   FILE *f;
   int more;
   int i;

   static int counter = 0; /* This supports animation sequences */
   char fname[64];


   //glViewport(0, 0, width, height);



   glUseProgram(program);
   glUniform1f(uniform_fade, cur_fade);

   printf("Generating %d by %d image file...\n", IMAGE_WIDTH, IMAGE_HEIGHT);

   /* allocate buffer large enough to store one tile */
   tile = malloc(TILE_WIDTH * TILE_HEIGHT * 3 * sizeof(GLubyte));
   if (!tile) {
      printf("Malloc of tile buffer failed!\n");
      return;
   }

   /* allocate buffer to hold a row of tiles */
   buffer = malloc(IMAGE_WIDTH * TILE_HEIGHT * 3 * sizeof(GLubyte));
   if (!buffer) {
      free(tile);
      printf("Malloc of tile row buffer failed!\n");
      return;
   }

   /* Setup.  Each tile is TILE_WIDTH x TILE_HEIGHT pixels. */
   tr = trNew();
   trTileSize(tr, TILE_WIDTH, TILE_HEIGHT, TILE_BORDER);
   trTileBuffer(tr, GL_RGB, GL_UNSIGNED_BYTE, tile);
   trImageSize(tr, IMAGE_WIDTH, IMAGE_HEIGHT);
   trRowOrder(tr, TR_TOP_TO_BOTTOM);

   trOrtho(tr, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

   // png output
   bit_depth = 8;
   // Open the file
   sprintf(fname,"supersize%04d.png",counter);
   FILE *fp = fopen(fname, "wb");
   if (!fp)
        abort_("[write_png_file] File %s could not be opened for writing", fname);
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
   info_ptr = png_create_info_struct(png_ptr);
   setjmp(png_jmpbuf(png_ptr));
   png_init_io(png_ptr, fp);
   setjmp(png_jmpbuf(png_ptr));
   png_set_IHDR(png_ptr, info_ptr, IMAGE_WIDTH, IMAGE_HEIGHT,
                 bit_depth, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        /* one metre = 100 centimetre (cm), and 2.54 cm = 1 inch */
		/* 1 metre is about 40 inches (well, 100/2.54 or 39.37) */
		/* so the number of dots per metre is about 40 times */
		/* larger than the number of dots per inch */
		/* thus DPM = DPI * 100 / 2.54 = DPI * 10000 / 254 */
		int ppm_x, ppm_y; /* pixels per metre */
        int dpi = 300;
		ppm_x = (dpi * 10000 + 127) / 254; /* round to nearest */
		ppm_y = ppm_x;
		png_set_pHYs(png_ptr, info_ptr, ppm_x, ppm_y,
			PNG_RESOLUTION_METER);

    png_write_info(png_ptr, info_ptr);
    setjmp(png_jmpbuf(png_ptr));


   /*
    * Should set GL_PACK_ALIGNMENT to 1 if the image width is not
    * a multiple of 4, but that seems to cause a bug with some NVIDIA
    * cards/drivers.
    */
   glPixelStorei(GL_PACK_ALIGNMENT, 1); // was 4

   /* Draw tiles */
   more = 1;
   while (more) {
      int curColumn;
      trBeginTile(tr);
      curColumn = trGet(tr, TR_CURRENT_COLUMN);
      draw_tr();      /* draw our stuff here */
      more = trEndTile(tr);

      /* save tile into tile row buffer*/
      {
	 int curTileWidth = trGet(tr, TR_CURRENT_TILE_WIDTH);
	 int bytesPerImageRow = IMAGE_WIDTH*3*sizeof(GLubyte);
	 int bytesPerTileRow = (TILE_WIDTH-2*TILE_BORDER) * 3*sizeof(GLubyte);
	 int xOffset = curColumn * bytesPerTileRow;
	 int bytesPerCurrentTileRow = (curTileWidth-2*TILE_BORDER)*3*sizeof(GLubyte);
	 int i;
	 int curTileHeight = trGet(tr, TR_CURRENT_TILE_HEIGHT);
	 for (i=0;i<curTileHeight;i++) {
	    memcpy(buffer + i*bytesPerImageRow + xOffset, /* Dest */
		   tile + i*bytesPerTileRow,              /* Src */
		   bytesPerCurrentTileRow);               /* Byte count*/
	 }
      }

      if (curColumn == trGet(tr, TR_COLUMNS)-1) {
	 /* write this buffered row of tiles to the file */
	 int curTileHeight = trGet(tr, TR_CURRENT_TILE_HEIGHT);
	 int bytesPerImageRow = IMAGE_WIDTH*3*sizeof(GLubyte);
	 int i;
	 GLubyte *rowPtr;
         /* The arithmetic is a bit tricky here because of borders and
          * the up/down flip.  Thanks to Marcel Lancelle for fixing it.
          */
	 for (i=2*TILE_BORDER;i<curTileHeight;i++) {
	    /* Remember, OpenGL images are bottom to top.  Have to reverse. */
	    rowPtr = buffer + (curTileHeight-1-i) * bytesPerImageRow;
        png_write_row(png_ptr, rowPtr);
	 }
      }

   }
   trDelete(tr);

   setjmp(png_jmpbuf(png_ptr));
   png_write_end(png_ptr, NULL);

   fclose(fp);
   printf("%s complete.\n", fname);
   counter++;

   free(tile);
   free(buffer);


          // put out view back to how it was originally
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
        glOrtho(-1.0, 1.0, -1.0, 1.0, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);

        // update it back to the original
        onDisplay();

   //exit(0);
}





void free_resources()
{
  glDeleteProgram(program);
  glDeleteBuffers(1, &vbo_triangle);
  //glDeleteBuffers(1, &vbo_triangle_colors);
}

void keys( unsigned char key, int x, int y )
{


    if( key == '1' ) {
        glClearColor(1.0, 1.0, 1.0, 1.0);
		//display();
        //printf("done init resources()\n");
	}
    if( key == '2' ) {
        glClearColor(0.0, 0.0, 0.0, 1.0);
		//update_resources();
        //printf("done init resources()\n");
	}

    if( key == '3' ) {
        glClearColor(0.1, 0.1, 0.1, 1.0);
		//update_resources();
        //printf("done init resources()\n");
	}




	if( key == ' ' ) {
		update_resources();
        //printf("done init resources()\n");
	}

    if( key == '='){
        nPts += 10;
        update_resources();
        //printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);
    }

    if( key == '-'){
        nPts -= 10;
        if(nPts<=40) nPts = 40;
        update_resources();
        //printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);

    }

    if( key == 'F'){
        cur_fade += 0.01;
        if(cur_fade > 1.0) cur_fade = 1.0;
    }

    if( key == 'f'){
        cur_fade -= 0.01;
        if(cur_fade < 0.0) cur_fade = 0.0;
    }

    if( key == ']'){
        edge_precent += 0.01;
        if(edge_precent > 1.0) edge_precent = 1.0;
        update_resources();
        //printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);
    }

    if( key == '['){
        edge_precent -= 0.01;
        if(edge_precent < 0.0) edge_precent = 0.0;
        update_resources();
        //printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);

    }

    if( key == 'p'){
        min_edge_value += 1.0;
        if(min_edge_value >= 254) min_edge_value = 254;
        update_resources();
        //printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);
    }

    if( key == 'o'){
        min_edge_value -= 1.0;
        if(min_edge_value < 0) min_edge_value = 0;
        update_resources();
        //printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);

    }

    if( key == 's'){
        WindowDump_PNG();
        //printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);

    }

    if( key == 'S'){
        TILE_WIDTH = glutGet(GLUT_WINDOW_WIDTH);
        TILE_HEIGHT = glutGet(GLUT_WINDOW_HEIGHT);
        tile_Display();
        //printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);

    }

    if( key == 'b'){
        do_border = !do_border;
        update_resources();
        //printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);

    }

    if( key == 'r'){
        do_random = !do_random;
        update_resources();
        //printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);

    }

    if( key == ','){
        area_threshold -= 100.0;
        if(area_threshold < 0) area_threshold = 0.0;

        update_resources();
    }

    if( key == '.'){
        area_threshold += 100.0;

        update_resources();
    }

    if(key == 't'){
        do_wireframe = !do_wireframe;
        //update_resources();
        update_resources_wireframe();
    }

    if(key == 'm'){
        do_textureMap = !do_textureMap;
        update_resources();
        //update_resources_wireframe();
    }

    if(key == 'g'){
        do_grayScale = !do_grayScale;
        update_resources();
        //update_resources_wireframe();
    }

    //If the user presses q
    if( key == 'a' )
    {
		//Cycle alias mode
		switch( gMode )
		{
			case ALIAS_MODE_ALIASED:
				printf( "Antialiased\n" );
				gMode = ALIAS_MODE_ANTIALIASED;
				break;

			case ALIAS_MODE_ANTIALIASED:
				printf( "Multisampled\n" );
				gMode = ALIAS_MODE_MULTISAMPLE;
				break;

			case ALIAS_MODE_MULTISAMPLE:
				printf( "Aliased\n" );
				gMode = ALIAS_MODE_ALIASED;
				break;
		}
    }

    printf("nPts = %d, cur_fade = %f, area_thresh = %f, edge_precent = %f, min_edge_value = %d\n", nPts, cur_fade, area_threshold, edge_precent, min_edge_value);
    onDisplay();
}




int main(int argc, char* argv[]) {

    GLint buf, sbuf;

  int     i,j,x,y;
  int   count;  // only for printf's debug
  int   got_color;

    float     *points;
    float   fac = RAND_MAX;
    int     nDim = 2;
    int     cwise = 1;
    pts     *point_list;

    unsigned int    *triIndexList;

    double  rnum;
    double  r1,r2;

    int  l, t;
    int nEdgePixels;

    png_byte* ptrp;
	png_byte* rowp;



    if (argc < 2){
        printf("usage: ./app input.png [output.png]\n");
        exit(1);
    }

    printf("nPts = %d, cur_fade = %f\n", nPts, cur_fade);


    // read and write a test png file
  read_png_file(argv[1]);
  // malloc memory for the edges
  edges = malloc2d_int(height,width);

  // set output supersized image size
  // originals
  IMAGE_WIDTH= 22 * width;
  IMAGE_HEIGHT= 22 * height;



  // add noise to image
  //add_noise();

  // get image edges
  sobel(edges);



  srand(time(NULL));
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA|GLUT_ALPHA|GLUT_DOUBLE|GLUT_DEPTH|GLUT_MULTISAMPLE);
  glutInitWindowSize(width, height);
  glutCreateWindow("TFU");

  GLenum glew_status = glewInit();
  if (glew_status != GLEW_OK) {
    fprintf(stderr, "Error: %s\n", glewGetErrorString(glew_status));
    return 1;
  }

    glGetIntegerv(GL_SAMPLE_BUFFERS, &buf);
    printf("number of sample buffers is %d\n", buf);
    glGetIntegerv(GL_SAMPLES, &sbuf);
    printf("number of samples is %d\n", sbuf);

  if (!GLEW_VERSION_2_0) {
    fprintf(stderr, "Error: your graphic card does not support OpenGL 2.0\n");
    return 1;
  }

    points = malloc(nDim*nPts*sizeof(float));
    point_list = malloc(nPts * sizeof(pts));

    l = 0;

    for(i=0;i<=10;i+=2){
        points[l] = (float)i/10.0; l++;
        points[l] = 0.0; l++;
    }

    for(i=1;i<=10;i+=2){
        points[l] = (float)i/10.0; l++;
        points[l] = 1.0; l++;
    }

    for(i=0;i<=10;i+=2){
        points[l] = 0.0; l++;
        points[l] = (float)i/10.0; l++;
    }

    for(i=0;i<=10;i+=2){
        points[l] = 1.0; l++;
        points[l] = (float)i/10.0; l++;
    }

    // how may pixels are in the image that represent an edge?
    //l = 0;
    do{
        r1 = rand_in_range(0.0,1.0);
        r2 = rand_in_range(0.0,1.0);
        // convert random numbers to image coordinates
        y = floor((float)height * r1);
        x = floor((float)width  * r2);
        //printf("testing pixel [%d][%d]: image pixel is %d\n",y,x,edges[y][x]);
        // is this on an adge?
        if(edges[y][x] > min_edge_value){   // found an edge
            //printf("\tfound an edge at [%d][%d]\n",y,x);
            points[l] = r2; l++;
            points[l] = r1; l++;
            // add edge pixel to out image for testing
            rowp = out_pointers[y];
            ptrp = &(rowp[(x)*4]);

            ptrp[0] = 255;
            ptrp[1] = 0;
            ptrp[2] = 0;
            ptrp[3] = 255;

        }

    }while(l<(nDim*nPts)*edge_precent);

    for(i=l;i<2*nPts;i++){
        points[i] = rand_in_range(0.0,1.0);
    }

    write_png_file(argv[2]);

    // need to read the original image again because we screwed it up in the sobol function
  // can fix this by not messing with the input png values (do the grayscale conversion to another 2d array)
  read_png_file(argv[1]);

    /*
    points[0] = 0; points[1] = 0;
    points[2] = 0; points[3] = 1;
    points[4] = 1; points[5] = 0;
    points[6] = 1; points[7] = 1;
    */






    l = 0;
    for(i=0;i<2*nPts;i+=2){
        point_list[l].x = points[i];
        point_list[l].y = points[i+1];
        l++;
    }


    //printf("#building triangle list\n");
    triIndexList = BuildTriangleIndexList((void*)points, fac, nPts, nDim, cwise, &nTriVert);
    //printf("#done\n");

    //printf("#numTriangle Vert = %d\n", nTriVert);
    //printf("#nTriangles = %d\n", nTriVert/3);


    /*
    l = 0;
    for(i=0;i<2*nPts;i+=2){
        printf("#point[ %d ]= %f,%f\n", l, points[i], points[i+1]);
        l++;
    }

    l = 0;
    for(i=0;i<nPts;i++){
        printf("#point_list[ %d ] = %f,%f\n", i, point_list[i].x, point_list[i].y);
        l++;
    }
    */



    // set up triangles for display
    nTri = (nTriVert/3);
    triangle_attributes = malloc(nTriVert * sizeof(attributes));
    count = 0;
    for(i=0;i<nTriVert;i+=3){
        //printf("#triangle %d hs vertices = %d, %d, %d\n",l, triIndexList[i],triIndexList[i+1],triIndexList[i+2]);
        /*
        r1 = rand_in_range(0.0, 0.2);
        r2 = rand_in_range(0.4, 1.0);
        rnum = rand_in_range(0.4, 1.0);

        triangle_attributes[i].v_color[0] = r1; //0.31640625 + r1;
        triangle_attributes[i].v_color[1] = r2; //0.65625 + r2;
        triangle_attributes[i].v_color[2] = rnum;

        triangle_attributes[i+1].v_color[0] = r1; //0.31640625 + r1;
        triangle_attributes[i+1].v_color[1] = r2; //0.65625 + r2;
        triangle_attributes[i+1].v_color[2] = rnum;

        triangle_attributes[i+2].v_color[0] = r1; //0.31640625 + r1;
        triangle_attributes[i+2].v_color[1] = r2; //0.65625 + r2;
        triangle_attributes[i+2].v_color[2] = rnum;
        */

        triangle_attributes[i].coord2d[0] = (point_list[triIndexList[i]].x - 0.5)/0.5;
        triangle_attributes[i].coord2d[1] = (point_list[triIndexList[i]].y - 0.5)/0.5;

        triangle_attributes[i+1].coord2d[0] = (point_list[triIndexList[i+1]].x - 0.5)/0.5;
        triangle_attributes[i+1].coord2d[1] = (point_list[triIndexList[i+1]].y - 0.5)/0.5;

        triangle_attributes[i+2].coord2d[0] = (point_list[triIndexList[i+2]].x - 0.5)/0.5;
        triangle_attributes[i+2].coord2d[1] = (point_list[triIndexList[i+2]].y - 0.5)/0.5;



        /*
        printf("%f %f\n", point_list[triIndexList[i]].x, point_list[triIndexList[i]].y);
        printf("%f %f\n", point_list[triIndexList[i+1]].x, point_list[triIndexList[i+1]].y);
        printf("%f %f\n", point_list[triIndexList[i+2]].x, point_list[triIndexList[i+2]].y);
        printf("%f %f\n", point_list[triIndexList[i]].x, point_list[triIndexList[i]].y);
        */

        /*
        printf("triangle %d:\n", count);
        printf("\t%f %f\n", triangle_attributes[i].coord2d[0], triangle_attributes[i].coord2d[1]);
        printf("\t%f %f\n", triangle_attributes[i+1].coord2d[0], triangle_attributes[i+1].coord2d[1]);
        printf("\t%f %f\n", triangle_attributes[i+2].coord2d[0], triangle_attributes[i+2].coord2d[1]);
        //printf("%f %f\n", triangle_attributes[i].coord2d[0], triangle_attributes[i].coord2d[1]);
        printf("\n");
        */

        count++;

    }
    //fflush(stdout);


    // process png file to get color for triangle vertex

    polygon = malloc(3*sizeof(pts));
    vertex_color = malloc(sizeof(vcolor));
    printf("image details: width = %d, height = %d\n", width, height);
    for(i=0;i<nTriVert;i+=3){   //for each triangle

        // rescale triangle coords to match image dimensions
        // points are in -1 <= x <= 1,     -1 <= y <= 1
        // image is       0 <= x <= width,  0 <= y <= height


        polygon[0].x = width * (0.5*triangle_attributes[i].coord2d[0] + 0.5);
        polygon[0].y = height *(0.5*triangle_attributes[i].coord2d[1] + 0.5);

        polygon[1].x = width * (0.5*triangle_attributes[i+1].coord2d[0] + 0.5);
        polygon[1].y = height *(0.5*triangle_attributes[i+1].coord2d[1] + 0.5);

        polygon[2].x = width * (0.5*triangle_attributes[i+2].coord2d[0] + 0.5);
        polygon[2].y = height * (0.5*triangle_attributes[i+2].coord2d[1] + 0.5);

        /*
        printf("polygon in image space:\n");
        printf("tri = %d\n", i/3);
        printf("\tvertex: %f, %f (was %f, %f)\n", polygon[0].x, polygon[0].y, 0.5*triangle_attributes[i].coord2d[0] + 0.5, 0.5*triangle_attributes[i].coord2d[1] + 0.5);
        printf("\tvertex: %f, %f (was %f, %f)\n", polygon[1].x, polygon[1].y, 0.5*triangle_attributes[i+1].coord2d[0]+0.5,0.5*triangle_attributes[i+1].coord2d[1]+0.5);
        printf("\tvertex: %f, %f (was %f, %f)\n", polygon[2].x, polygon[2].y, 0.5*triangle_attributes[i+2].coord2d[0]+0.5, 0.5*triangle_attributes[i+2].coord2d[1]+0.5);
        */

        got_color = process_triangle_poly(polygon, vertex_color);    // find out what pixels are in this triangle

        if(got_color == 1){
            triangle_attributes[i].v_color[0] = vertex_color->r;
            triangle_attributes[i].v_color[1] = vertex_color->g;
            triangle_attributes[i].v_color[2] = vertex_color->b;

            triangle_attributes[i+1].v_color[0] = vertex_color->r;
            triangle_attributes[i+1].v_color[1] = vertex_color->g;
            triangle_attributes[i+1].v_color[2] = vertex_color->b;

            triangle_attributes[i+2].v_color[0] = vertex_color->r;
            triangle_attributes[i+2].v_color[1] = vertex_color->g;
            triangle_attributes[i+2].v_color[2] = vertex_color->b;
        }
        else{   // actually want the color from a neighbor triangle
                // i don't think we have an adjacency list. damn!
            triangle_attributes[i].v_color[0] = triangle_attributes[i-3].v_color[0];
            triangle_attributes[i].v_color[1] = triangle_attributes[i-3].v_color[0];
            triangle_attributes[i].v_color[2] = triangle_attributes[i-3].v_color[0];

            triangle_attributes[i+1].v_color[0] = triangle_attributes[i-2].v_color[0];
            triangle_attributes[i+1].v_color[1] = triangle_attributes[i-2].v_color[0];
            triangle_attributes[i+1].v_color[2] = triangle_attributes[i-2].v_color[0];

            triangle_attributes[i+2].v_color[0] = triangle_attributes[i-1].v_color[0];
            triangle_attributes[i+2].v_color[1] = triangle_attributes[i-1].v_color[0];
            triangle_attributes[i+2].v_color[2] = triangle_attributes[i-1].v_color[0];

        }
    }


    //write_png_file(argv[2]);



  if (init_resources()) {
    glutDisplayFunc(onDisplay);
    //glutDisplayFunc(tile_Display);
    glutIdleFunc(onIdle);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    glDisable( GL_LINE_SMOOTH );
    glDisable( GL_POLYGON_SMOOTH );
    glEnable( GL_MULTISAMPLE );

    glutKeyboardFunc( keys );
    glutMainLoop();
  }

  free_resources();
  return 0;
}
